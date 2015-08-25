#include "BBSMPSTree.hpp"

using namespace std;

double intTol=1e-6;


// Outputs solver status:
void outputLPStatus(solverState lpStatus) {
	
	string status;
	switch(lpStatus) {
		case Uninitialized:
		status = "Uninitialized";
		break;
		case LoadedFromFile:
		status = "LoadedFromFile";
		break;
		case Initialized:
		status = "Initialized";
		break;
		case PrimalFeasible:
		status = "PrimalFeasible";
		break;
		case DualFeasible:
		status = "DualFeasible";
		break;
		case Optimal:
		status = "Optimal";
		break;
		case ProvenUnbounded:
		status = "ProvenUnbounded";
		break;
		case ProvenInfeasible:
		status = "ProvenInfeasible";
		break;
		case Stopped:
		status = "Stopped";
		break;
	}
	BBSMPS_ALG_LOG_SEV(info) <<"PIPS-S has returned "<< status;
}



// Overload the "less than" operator so that priority_queue can use it
// for heapifying comparisons. Make BBSMPSNode a templated
// class when adding additional branching heuristics?
bool operator< (const BBSMPSNode& left,
	const BBSMPSNode&right)
{
	return left.getParentObjective() < right.getParentObjective();// left.parentObj < right.parentObj;
}

bool operator> (const BBSMPSNode& left,
	const BBSMPSNode& right)
{
	return  left.getParentObjective() > right.getParentObjective();//left.parentObj > right.parentObj;
}


BBSMPSTree::BBSMPSTree(const SMPSInput& smps): 
objUB(COIN_DBL_MAX),
objLB(-COIN_DBL_MAX),
optGapTol(1e-6),
lpPrimalTol(1e-6),
lpDualTol(1e-6),
compTol(lpPrimalTol),
status(LoadedFromFile),
nodesel(BestBound)
{

	BBSMPSSolver::initialize(smps); 

    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Calling B&B tree constructor.";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start

    
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    int mype=BBSMPSSolver::instance()->getMype();

    rootSolver.setPrimalTolerance(lpPrimalTol);
    rootSolver.setDualTolerance(lpDualTol);

    // TODO: Replace with real presolve.
    // For now, "cheat" by solving root LP before populating root node.
    // A warm start of the root node means that the root node solve
    // inside the B&B tree does not cost very much -- just the PIPS-S overhead.
    // However, this step is necessary in order to properly allocate primal
    // variables (due to slacks) and to determine which variables are basic.
    // This step is also currently required to instantiate lower & upper
    // bounds. In theory, the bounds could be obtained from the SMPS file.
    // In practice, getting the bounds from the SMPS file is cumbersome,
    // because first stage bounds and second stage bounds must be
    // queried separately and are returned as std::vector<double>s.
    // In addition, distributed data structures dictate some care in
    // how the assignments are performed: there must be checks to
    // ensure that the data to be assigned is owned by the "right" process.

    rootSolver.go();

    // Get lower & upper bounds on decision variables in LP.
    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Getting bounds for root node from presolve.";
    denseBAVector lb(rootSolver.getLB()), ub(rootSolver.getUB());

       // Allocate current best primal solution; normally this primal solution
    // is for the upper bound, but here, we have only the solution to an
    // LP relaxation, which may not be primal feasible. We don't check
    // primal/integer feasibility here.
    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Allocating primal solution.";
    ubPrimalSolution.allocate(dimsSlacks, ctx, PrimalVector);
    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Getting primal solution";
    ubPrimalSolution.copyFrom(rootSolver.getPrimalSolution());
    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "MIP Primal solution updated"; 

    // Update state of primal variables + slacks; slacks must be included
    // because these are used in a reformulation of the problem to standard
    // form: Ax + s = b, s >= 0.
    BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
    rootSolver.getStates(states);

    // Update global lower bound; really need to check feasibility, etc.
    // here, but I'm going to move this code to the B&B tree.
    double lpObj = rootSolver.getObjective();
    if ((lpObj - compTol) >= objLB) objLB = lpObj;

    BBSMPSNode *rootNode= new BBSMPSNode(lpObj,states);

    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Pushing root node onto B&B tree.";
    heap.push(rootNode);
    //if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Exiting B&B constructor.";

    BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
    branchingRuleManager.addBranchingRule(mfbr);


    BBSMPSHeuristicRounding *hr= new BBSMPSHeuristicRounding(1,1,"SimpleRounding");
    heuristicsManager.addHeuristic(hr);
    bbIterationCounter=0;
}



int BBSMPSTree::getFirstStageMinIntInfeasCol(const denseBAVector& primalSoln) {
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	int col;

		    // Return first index of integer variable with fractional value
	for (col = 0; col < input.nFirstStageVars(); col++)
	{

		bool isColInteger = input.isFirstStageColInteger(col);
		bool isValInteger = isIntFeas(primalSoln.getFirstStageVec()[col],
			intTol);

			// If the (col)th 1st stage primal variable is integer,
			// but has a fractional value, return idx
		if(isColInteger && !isValInteger) return col;
	}

		    // Otherwise, 1st stage is integer feasible: return -1;
	return -1;
}




		// NOTE: MPI standard requires passing ints, not bools
int BBSMPSTree::isFirstStageIntFeas(const denseBAVector& primalSoln) {
	return (getFirstStageMinIntInfeasCol(primalSoln) == -1);
}

int BBSMPSTree::getSecondStageMinIntInfeasCol(const denseBAVector& primalSoln, int scen) {
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    // Only makes sense when called by process that owns scenario scen
	assert(ctx.assignedScenario(scen));

	int col;
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col],
			intTol);

			// If the (col)th 2nd stage primal variable of the (scen)th
			// scenario is integer, but has fractional value, return idx
		if (isColInteger && !isValInteger) return col;
	}

		    // Otherwise, return -1;
	return -1;
}

int BBSMPSTree::isSecondStageIntFeas(const denseBAVector& primalSoln, int scen) {
	return (getSecondStageMinIntInfeasCol(primalSoln, scen) == -1);
}

bool BBSMPSTree::isLPIntFeas(const denseBAVector& primalSoln) {

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int is1stStageIntFeas(isFirstStageIntFeas(primalSoln));

		    // Check to see if all 2nd stage scenarios on current MPI rank are
		    // integer feasible. At first scenario that is integer infeasible,
		    // we know that the primal solution is not integer feasible.
	int isMyRankIntFeas(1);
	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			if(!isSecondStageIntFeas(primalSoln, scen)) {
				isMyRankIntFeas = 0;
				break;
			}
		}
	}

	int is2ndStageIntFeas(0);
	int errorFlag = MPI_Allreduce(&isMyRankIntFeas,
		&is2ndStageIntFeas,
		1,
		MPI_INT,
						  MPI_LAND, // MPI logical and
						  ctx.comm());
		    // TODO: Some error handling here

	return (is1stStageIntFeas && is2ndStageIntFeas);
}




		  // TODO: Add way to access single scenario.

		  /* TODO: Refactor MIP solver status into its own class to handle
		     the transition logic. */
		  /* TODO: Add more nuanced version of solver statuses. */
		  /* Methods defining allowable transitions for status */
		  /* status transition state diagram */
		  // NOTE: "Bounded" is not currently a solver state;
		  // TODO: Add Bounded as a solver state.
		  // LoadedFromFile -> {Bounded, PrimalFeasible, Optimal,
		  //                    ProvenInfeasible, Unbounded, Stopped}
		  // Bounded -> {PrimalFeasible, Optimal, ProvenInfeasible, Stopped}
		  // PrimalFeasible -> {Optimal, Stopped}

		  // void setStatusToBounded() {
		  //   bool isInReachableState = (LoadedFromFile == status);
		  //   if (isInReachableState) {
		  //     status = Bounded;
		  //   }
		  // }



void BBSMPSTree::branchAndBound() {
	int mype=BBSMPSSolver::instance()->getMype();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();

	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Starting branch-and-bound.";

	/* While heap not empty and there are still nodes in tree */
	// TODO: Add tolerance on optimality gap, time limit option.
	while (true) {

		/* If heap is empty, update status to Stopped (if possible) and break. */
		if (heap.empty()) {
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
			status.setStatusToStopped();

		/* If solver status is primal feasible, and the heap is empty, then
		the solution must be optimal. */
			if (status.isPrimalFeasible()) {
				status.setStatusToOptimal();
				if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Optimal solution found.";
			}

		/* If solver status is not primal feasible, then the MILP must be
		infeasible. */
		// TODO: Add test for unboundedness.
			if (status.isLoadedFromFile()) {
				status.setStatusToProvenInfeasible();
				if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "MILP is infeasible.";
			}

			break;
		}

		/* Get top-most node and pop it off of heap. */
		BBSMPSNode *currentNode_ptr=(heap.top());
		if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Copying node " << currentNode_ptr->getNodeNumber() << " off tree.";

		//if (0 == mype) cout << "Popping node " << nodeNumber << " off tree!\n";
		heap.pop();


		/* Set bounds of LP decision variables from BBSMPSNode */
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Setting bounds for LP subproblem.";
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Parent objective of this node "<< currentNode_ptr->getParentObjective();
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Parent pointer of this node "<< (currentNode_ptr->getParentPtr()!=NULL);

		denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
		denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());

		currentNode_ptr->getAllBranchingInformation(lb,ub);

		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		/* Set information on basic/nonbasic variables for warm starting */
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Setting warm start information.";
		BAFlagVector<variableState> ps;
		currentNode_ptr->getWarmStartState(ps);

		rootSolver.setStates(ps);
		rootSolver.commitStates();
		/* Solve LP defined by current node*/
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Solving LP subproblem.";
		rootSolver.go();

		/* Check solver status for infeasibility/optimality */
		solverState lpStatus = rootSolver.getStatus();

		// Only realistic solver states upon completion:
		// ProvenInfeasible, Optimal, ProvenUnbounded
		// Other solver states are intermediate states that should not
		// hold upon return from rootSolver.
		if (0 == mype) outputLPStatus(lpStatus);

		bool isLPinfeasible = (ProvenInfeasible == lpStatus); 
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "isLPinfeasible = " << isLPinfeasible;
		bool isLPunbounded = (ProvenUnbounded == lpStatus);
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "isLPunbounded = " << isLPunbounded ;
		bool isLPoptimal = (Optimal == lpStatus);
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "isLPoptimal = " << isLPoptimal;
		bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "isLPother = " << isLPother;
		assert (!isLPother); // Error if not infeasible/unbounded/optimal

		/* Fathom by infeasibility */
		// If LP solver returns infeasibility, fathom node, go to start of loop
		if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Checking for infeasibility...";
		if (isLPinfeasible) {
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by infeasibility.";
			delete currentNode_ptr;
			continue;
		}

		// Otherwise, LP is feasible. LP may be optimal or unbounded.
		// If LP is unbounded, so is the MILP.
		if (isLPunbounded) {
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "LP relaxation of node " << currentNode_ptr->getNodeNumber()
				<< " is unbounded.\n"
			<< "Please add additional constraints to "
			<< "bound the MILP.";
			return;
		}

		// At this point, LP must be optimal.
		assert (isLPoptimal); // Error if not optimal.

		// TODO: Combine the integrality and branching steps later

		/* If LP solution is optimal, can fathom by value dominance. */

		if (isLPoptimal) {
			// If LP solver returns optimal, then the objective is bounded below.
			// TODO: Change solver status to "Bounded".
			//	setStatusToBounded();

			/* Lower bound update */ // This part is still incorrect!
			// If the LP solver returns an optimal solution AND that solution is
			// greater than the current best lower bound, update the best lower bound
			// on the objective function value.

			// Since the branch-and-bound tree is stored as a min-heap, the current
			// node being explored always has the minimal objective function value.
			// If the branch-and-bound tree is ever re-heapified so that it is NOT
			// a min-heap, but has the heap property for some other ordering, then
			// this update cannot occur without taking the min objective function
			// value over all values of the parent objective function for nodes
			// still in the B&B tree (which is expensive if it is not the key used
			// to heapify the min-heap).
			//objLB = currentNode.parentObj;
			//if ((lpObj - compTol) >= objLB) {
			//  if (0 == mype) cout << "Current best lower bound is " << objLB << endl;
			//  if (0 == mype) cout << "Updating best lower bound to " << lpObj << endl;
			//  objLB = lpObj;
			//}

			/* Fathom by value dominance */
			// Optimal LP objective function value is lower bound on the objective
			// function value of the LP derived from any node in the subtree
			// of the B&B tree rooted at the current node, so no feasible
			// solution in that subtree can have a lesser objective function
			// value than the current upper bound on the optimal value of
			// the MILP objective function.

			// Get LP objective function value.
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Getting LP objective...";
			double lpObj = rootSolver.getObjective();

			currentNode_ptr->setObjective(lpObj);

			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Checking for value dominance...";

			if ((lpObj - compTol) >= objUB) {
				if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by value dominance.";

				delete currentNode_ptr;

				if (nodesel == BestBound) {
					if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Can stop if best bound node selection rule";
					status.setStatusToOptimal();
					if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "All nodes can be fathomed! Terminating.";
					break;
				}
				continue;
			}
		}

		/* Get primal solution */
		//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Getting primal solution...";
		denseBAVector primalSoln(rootSolver.getPrimalSolution());


		vector<denseBAVector> heuristicSolutions;
		heuristicsManager.runHeuristics(currentNode_ptr,primalSoln,heuristicSolutions);
		if (heuristicSolutions.size()>0) {
			bool isFeas=isLPIntFeas(heuristicSolutions[0]);
			BBSMPS_ALG_LOG_SEV(info)<<"Heuristic found solution.";
		}
			/* If primal solution is integral: */
			//  - Update solver status to PrimalFeasible
			//  - Check if upper bound improved
			//  - If so, update current primal solution for that upper bound
			//  - Fathom node, go to start of loop

		if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Checking for integrality of primal solution...";
		if(isLPIntFeas(primalSoln)) {
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Node " << currentNode_ptr->getNodeNumber() << " is integer feasible.";
			status.setStatusToPrimalFeasible();

			// TODO: Maintain solution pool of best k solutions for k = some small value

			/* Update upper bound if it's less than current best upper bound, and
			the LP solution is optimal (not unbounded). */
			double newUB = rootSolver.getObjective();
			bool isNewUBbetter = (newUB < (objUB - compTol));
			if (isLPoptimal && isNewUBbetter) {
				if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Updating best upper bound to " << newUB ;
				objUB = rootSolver.getObjective();
				ubPrimalSolution.copyFrom(primalSoln);
			}
			delete currentNode_ptr;
			continue;
		}

		// TODO: Fathom by value dominance in breadth-first fashion?

		/* Otherwise, primal solution is not integral: */
		// - Check to see if lower bound can be updated
		// - Branch

		/* Optimality gap termination criteria: if gap between best
		upper bound and best lower bound on objective function is
		less than the optimality gap, stop the solver and return the
		feasible solution corresponding to the best upper bound. */
		if (abs(objUB - objLB) <= optGapTol) {
			assert(objLB <= objUB); // this statement could be tripped by numerical error
			status.setStatusToOptimal();
			if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Optimality gap reached! Terminating.";
			break;
		}

		double lpObj = rootSolver.getObjective();
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
		BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
		BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
		rootSolver.getStates(states);

			//currentNode_ptr->setWarmStartState(states);


		vector<BBSMPSNode*> children;

		branchingRuleManager.branch(currentNode_ptr,children,primalSoln);

		if (children.size()>0){
			for (int i=0; i<children.size();i++){
				children[i]->setWarmStartState(states);
				heap.push(children[i]);
			}

		}
		bbIterationCounter++;
		if (0 == mype) {
			BBSMPS_ALG_LOG_SEV(info)<<"\n----------------------------------------------------\n"<<
			"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":Tree Size:"<<heap.size()<<"\n"<<
			"----------------------------------------------------";
		}
	}

//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Objective function value = " << objUB ;
//if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Objective function LB = " << objLB ;
	if (0 == mype) {
		BBSMPS_ALG_LOG_SEV(info)<<"\n--------------EXPLORATION TERMINATED----------------\n"<<
		"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":Tree Size:"<<heap.size()<<"\n"<<
		"----------------------------------------------------";
		heuristicsManager.printStatistics();
		branchingRuleManager.printStatistics();
	}

}
