/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
#include "BBSMPSTree.hpp"

using namespace std;

double intTol=1e-6;
BBSMPSNode *BBSMPSTree::rootNode =NULL;

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
nodesel(BestBound),
tiLim(COIN_INT_MAX),
nodeLim(COIN_INT_MAX),
solsDiscoveredLimit(COIN_INT_MAX),
solsDiscoveredInit(0),
verbosityActivated(true),
cuttingPlanesManager(5)
{
	double timeStart=MPI_Wtime();

	BBSMPSSolver::initialize(smps);
	double timeStampPreProc=MPI_Wtime();
	PreProcessingTime=timeStampPreProc-timeStart;
   // if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Calling B&B tree constructor.";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start



    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
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

   // BBSMPSCuttingPlane bcp;
		//for (int i=0; i<=bbIterationCounter; i++)
	//bcp.applyCuttingPlane();

    LPRelaxationTime=MPI_Wtime()-timeStampPreProc;
    // Get lower & upper bounds on decision variables in LP.
    if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting bounds for root node from presolve.";
    denseBAVector lb(rootSolver.getLB()), ub(rootSolver.getUB());
       // Allocate current best primal solution; normally this primal solution
    // is for the upper bound, but here, we have only the solution to an
    // LP relaxation, which may not be primal feasible. We don't check
    // primal/integer feasibility here.
    //pwdif (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Allocating primal solution.";
    ubPrimalSolution.allocate(dimsSlacks, ctx, PrimalVector);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting primal solution";
    ubPrimalSolution.copyFrom(rootSolver.getPrimalSolution());
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "MIP Primal solution updated";
    BBSMPSSolver::instance()->setLPRelaxation(ubPrimalSolution);
    // Update state of primal variables + slacks; slacks must be included
    // because these are used in a reformulation of the problem to standard
    // form: Ax + s = b, s >= 0.

    BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);

    rootSolver.getStates(states);
    BBSMPSSolver::instance()->setOriginalWarmStart(states);
    // Update global lower bound; really need to check feasibility, etc.
    // here, but I'm going to move this code to the B&B tree.
    double lpObj = rootSolver.getObjective();
    LPRelaxationValue=lpObj;
    if ((lpObj - compTol) >= objLB) objLB = lpObj;

    std::vector< std::pair < BAIndex, variableState > > emptyStates;

    BBSMPSNode *rootNode= new BBSMPSNode(lpObj,emptyStates);

    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Pushing root node onto B&B tree.";
    heap.push(rootNode);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Exiting B&B constructor.";

    BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
   branchingRuleManager.addBranchingRule(mfbr);


    //BBSMPSPseudoCostBranchingRule *mfbr2= new BBSMPSPseudoCostBranchingRule(100);
  	//  branchingRuleManager.addBranchingRule(mfbr2);

    bbIterationCounter=0;
    nodesFathomed=0;
    nodesBecameInteger=0;
    BBSMPSTree::rootNode=rootNode;
}



BBSMPSTree::BBSMPSTree(BBSMPSNode *node, double lb, double ub):
objUB(ub),
objLB(lb),
optGapTol(1e-6),
lpPrimalTol(1e-6),
lpDualTol(1e-6),
compTol(lpPrimalTol),
status(LoadedFromFile),
nodesel(BestBound),
tiLim(COIN_INT_MAX),
nodeLim(COIN_INT_MAX),
solsDiscoveredLimit(COIN_INT_MAX),
solsDiscoveredInit(0),
verbosityActivated(true),
cuttingPlanesManager(2)
{

	//This initialization assumes the solver class has already been intialized
	assert(BBSMPSSolver::isInitialized());

    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Calling B&B tree constructor.";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start

    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();

    double lpObj = -INFINITY;
    if ((lpObj - compTol) >= objLB) objLB = lpObj;


    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Pushing root node onto B&B tree.";
    heap.push(node);
    //if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Exiting B&B constructor.";

    BBSMPSMaxFracBranchingRule *mfbr= new BBSMPSMaxFracBranchingRule(10);
    branchingRuleManager.addBranchingRule(mfbr);


    bbIterationCounter=0;
     nodesFathomed=0;
   nodesBecameInteger=0;
   removeCuts();
	//node->getAllCuttingUids(currentlyAppliedPlanes);



}

BBSMPSTree::~BBSMPSTree(){
	//Dismantling nodes
	while (!heap.empty()){
		BBSMPSNode *currentNode_ptr=(heap.top());
		currentNode_ptr->eliminate();
		heap.pop();

	}

	//Dismantling Branching Rules
  	branchingRuleManager.freeResources();

	//Dismantling Heuristics
	heuristicsManager.freeResources();

	//Dismantling Cutting Planes
	cuttingPlanesManager.freeResources();
}



void BBSMPSTree::generateIncrementalWarmState(BBSMPSNode* node, const BAFlagVector<variableState> & originalState, const BAFlagVector<variableState> &currentState){

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();

	std::vector< std::pair < BAIndex, variableState > > changes;


	const denseFlagVector<variableState> &stageVec = currentState.getFirstStageVec();
	const denseFlagVector<variableState> &stageVec2 = originalState.getFirstStageVec();
	for (int j = 0; j < stageVec.length(); j++) {
		if (stageVec[j]!=stageVec2[j]){
			BAIndex aux;
			aux.scen=-1;
			aux.idx=j;
			changes.push_back(std::pair < BAIndex, variableState > (aux,stageVec[j]));

		}
	}


	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			const denseFlagVector<variableState> &stageVec = currentState.getSecondStageVec(scen);
			const denseFlagVector<variableState> &stageVec2 = originalState.getSecondStageVec(scen);
			for (int j = 0; j < stageVec.length(); j++) {
				if (stageVec[j]!=stageVec2[j]){
					BAIndex aux;
					aux.scen=scen;
					aux.idx=j;
					changes.push_back(std::pair < BAIndex, variableState > (aux,stageVec[j]));

				}
			}
		}
	}



	node->setIncrementalWarmStartState(changes);



}

void BBSMPSTree::removeCuts(){


	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	const BADimensionsSlacks &originalDimsSlacks= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	rootSolver.setLB(lb);
	rootSolver.setUB(ub);

	bool modelChanged=false;

	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			if (originalDimsSlacks.inner.numSecondStageCons(scen) < dimsSlacks.inner.numSecondStageCons(scen)){
				rootSolver.deleteLastSecondStageConsecutiveRows(scen,dimsSlacks.inner.numSecondStageCons(scen)-originalDimsSlacks.inner.numSecondStageCons(scen));
				modelChanged=true;
			}
			if (originalDimsSlacks.inner.numSecondStageVars(scen) < dimsSlacks.inner.numSecondStageVars(scen)){
				rootSolver.deleteLastSecondStageConsecutiveColumns(scen,dimsSlacks.inner.numSecondStageVars(scen)-originalDimsSlacks.inner.numSecondStageVars(scen));
				modelChanged=true;
			}
		}

	}

	if (originalDimsSlacks.inner.numFirstStageCons() < dimsSlacks.inner.numFirstStageCons()){
		rootSolver.deleteLastFirstStageConsecutiveRows(dimsSlacks.inner.numFirstStageCons()-originalDimsSlacks.inner.numFirstStageCons());
		modelChanged=true;
	}
	if (originalDimsSlacks.inner.numFirstStageVars() < dimsSlacks.inner.numFirstStageVars()){
		rootSolver.deleteLastFirstStageConsecutiveColumns(dimsSlacks.inner.numFirstStageVars()-originalDimsSlacks.inner.numFirstStageVars());
		modelChanged=true;
	}
	if (modelChanged){
		BBSMPSSolver::instance()->commitNewColsAndRows();
	}
	currentlyAppliedPlanes.clear();

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



	/* While heap not empty and there are still nodes in tree */
	// TODO: Add tolerance on optimality gap, time limit option.
	while (true) {

		PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
		if ((BBSMPSSolver::instance()->getWallTime())>tiLim){
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Time Limit reached.";
			status.setStatusToStopped();
			break;
		}
		if (BBSMPSSolver::instance()->getSolPoolSize()-solsDiscoveredInit>=solsDiscoveredLimit){
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Solution Limit reached.";
			status.setStatusToStopped();
			break;
		}

		if (bbIterationCounter>nodeLim){
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node Limit reached.";
			status.setStatusToStopped();
			break;
		}
		/* If heap is empty, update status to Stopped (if possible) and break. */
		if (heap.empty()) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Heap is empty.";
			status.setStatusToStopped();

		/* If solver status is primal feasible, and the heap is empty, then
		the solution must be optimal. */
			if (status.isPrimalFeasible()) {
				status.setStatusToOptimal();
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimal solution found.";
			}

		/* If solver status is not primal feasible, then the MILP must be
		infeasible. */
		// TODO: Add test for unboundedness.
			if (status.isLoadedFromFile()) {
				status.setStatusToProvenInfeasible();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "MILP is infeasible.";
			}

			break;
		}

		/* Get top-most node and pop it off of heap. */
		BBSMPSNode *currentNode_ptr=(heap.top());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Copying node " << currentNode_ptr->getNodeNumber() << " off tree.";

		heap.pop();

		if (nodesel == BestBound) {
			objLB=currentNode_ptr->getParentObjective();
			if ((objLB - compTol) >= objUB) {
				objLB=objUB;
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Can stop if best bound node selection rule";
				status.setStatusToOptimal();
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "All nodes can be fathomed! Terminating.";
				currentNode_ptr->eliminate();
				break;
			}

		}

		vector<int> nodeCuttingPlaneUids;
		currentNode_ptr->getAllCuttingUids(nodeCuttingPlaneUids);
		if (nodeCuttingPlaneUids.size()!=currentlyAppliedPlanes.size() || !equal(nodeCuttingPlaneUids.begin(), nodeCuttingPlaneUids.begin() + nodeCuttingPlaneUids.size(), currentlyAppliedPlanes.begin())){
			removeCuts();
			std::vector<BBSMPSCuttingPlane*> cpVector;
			currentNode_ptr->getAllCuttingPlanes(cpVector);
			for (int i=0; i< cpVector.size(); i++){
				cpVector[i]->applyCuttingPlane();
			}
			BBSMPSSolver::instance()->commitNewColsAndRows();
			currentlyAppliedPlanes=nodeCuttingPlaneUids;
		}

		/* Set bounds of LP decision variables from BBSMPSNode */
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds for LP subproblem.";
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent objective of this node "<< currentNode_ptr->getParentObjective();
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Parent pointer of this node "<< (currentNode_ptr->getParentPtr()!=NULL);

		denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
		denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "getting branching info.";

		currentNode_ptr->getAllBranchingInformation(lb,ub);
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Setting bounds.";

		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		/* Set information on basic/nonbasic variables for warm starting */

		BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());

		currentNode_ptr->reconstructWarmStartState(ps);

		rootSolver.setStates(ps);

		rootSolver.commitStates();

		/* Solve LP defined by current node*/
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Solving LP subproblem.";
		rootSolver.go();

		/* Check solver status for infeasibility/optimality */
		solverState lpStatus = rootSolver.getStatus();

		// Only realistic solver states upon completion:
		// ProvenInfeasible, Optimal, ProvenUnbounded
		// Other solver states are intermediate states that should not
		// hold upon return from rootSolver.
		if (0 == mype) outputLPStatus(lpStatus);

		bool isLPinfeasible = (ProvenInfeasible == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPinfeasible = " << isLPinfeasible;
		bool isLPunbounded = (ProvenUnbounded == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPunbounded = " << isLPunbounded ;
		bool isLPoptimal = (Optimal == lpStatus);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPoptimal = " << isLPoptimal;
		bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "isLPother = " << isLPother;
		assert (!isLPother); // Error if not infeasible/unbounded/optimal

		/* Fathom by infeasibility */
		// If LP solver returns infeasibility, fathom node, go to start of loop
		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for infeasibility...";
		if (isLPinfeasible) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by infeasibility.";
			currentNode_ptr->eliminate();
			nodesFathomed++;


			continue;
		}

		// Otherwise, LP is feasible. LP may be optimal or unbounded.
		// If LP is unbounded, so is the MILP.
		if (isLPunbounded) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "LP relaxation of node " << currentNode_ptr->getNodeNumber()
				<< " is unbounded.\n"
			<< "Please add additional constraints to "
			<< "bound the MILP.";

			return;
		}

		// At this point, LP must be optimal.
		assert (isLPoptimal); // Error if not optimal.

		denseBAVector primalSoln(rootSolver.getPrimalSolution());

		bool newCuttingPlanes=cuttingPlanesManager.generateCuttingPlanes(currentNode_ptr,primalSoln);

		if(newCuttingPlanes){
			currentNode_ptr->getAllCuttingUids(currentlyAppliedPlanes);
		}
		lpStatus = rootSolver.getStatus();
		primalSoln=denseBAVector(rootSolver.getPrimalSolution());
		isLPoptimal = (Optimal == lpStatus);
		ps=BAFlagVector<variableState>(BBSMPSSolver::instance()->getOriginalWarmStart());
		currentNode_ptr->reconstructWarmStartState(ps);
		// TODO: Combine the integrality and branching steps later

		/* If LP solution is optimal, can fathom by value dominance. */
		assert (isLPoptimal); // Error if not optimal.
		if (isLPoptimal) {
			// If LP solver returns optimal, then the objective is bounded below.
			// TODO: Change solver status to "Bounded".
			//	setStatusToBounded();

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
			//  if (0 == mype) //cout << "Current best lower bound is " << objLB << endl;
			//  if (0 == mype) //cout << "Updating best lower bound to " << lpObj << endl;
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
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Getting LP objective...";
			double lpObj = rootSolver.getObjective();

			currentNode_ptr->setObjective(lpObj);

			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for value dominance...";

			if ((lpObj - compTol) >= objUB) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Fathoming node " << currentNode_ptr->getNodeNumber() << " by value dominance.";

				currentNode_ptr->eliminate();
				nodesFathomed++;

				continue;
			}



		}


		/* Get primal solution */
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Getting primal solution...";

		denseBAVector originalSpaceSolution;
		originalSpaceSolution.allocate(BBSMPSSolver::instance()->getOriginalBADimensionsSlacks(), BBSMPSSolver::instance()->getBAContext(), PrimalVector);
		originalSpaceSolution.copyAndShrinkToDims(primalSoln);

		/* If primal solution is integral: */
		//  - Update solver status to PrimalFeasible
		//  - Check if upper bound improved
		//  - If so, update current primal solution for that upper bound
		//  - Fathom node, go to start of loop

		if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Checking for integrality of primal solution...";
		if(isLPIntFeas(primalSoln)) {
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Node " << currentNode_ptr->getNodeNumber() << " is integer feasible.";
			status.setStatusToPrimalFeasible();

			/* Update upper bound if it's less than current best upper bound, and
			the LP solution is optimal (not unbounded). */
			double newUB = rootSolver.getObjective();
			bool isNewUBbetter = (newUB < (objUB - compTol));
			if (isLPoptimal && isNewUBbetter) {
				if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Updating best upper bound to " << newUB ;
				objUB = rootSolver.getObjective();
				ubPrimalSolution.copyFrom(originalSpaceSolution);
				BBSMPSSolution aux(originalSpaceSolution,newUB, BBSMPSSolver::instance()->getWallTime());
				BBSMPSSolver::instance()->addSolutionToPool(aux);
			}
			currentNode_ptr->eliminate();


			nodesBecameInteger++;
			continue;
		}

		/* Optimality gap termination criteria: if gap between best
		upper bound and best lower bound on objective function is
		less than the optimality gap, stop the solver and return the
		feasible solution corresponding to the best upper bound. */

		double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
		if (gap <= optGapTol) {
			assert(objLB <= objUB); // this statement could be tripped by numerical error
			status.setStatusToOptimal();
			if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(info) << "Optimality gap reached! Terminating.";
			break;
		}

		double lpObj = rootSolver.getObjective();
		const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
		BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
		BAFlagVector<variableState> states(dimsSlacks, ctx, PrimalVector);
		rootSolver.getStates(states);

		vector<BBSMPSNode*> children;

		branchingRuleManager.branch(currentNode_ptr,children,originalSpaceSolution);

		if (children.size()>0){
			for (int i=0; i<children.size();i++){
				generateIncrementalWarmState(children[i], ps, states);
				children[i]->setObjective(lpObj);

				heap.push(children[i]);
			}

		}
		//removeCuts();
		heuristicsManager.runLPHeuristics(currentNode_ptr,primalSoln);


		if(heuristicsManager.willMIPHeuristicsRun(currentNode_ptr,originalSpaceSolution)){
			removeCuts();
			heuristicsManager.runMIPHeuristics(currentNode_ptr,originalSpaceSolution);
			//removeCuts();
		}



		if (BBSMPSSolver::instance()->getSolPoolSize()>0) objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

		bbIterationCounter++;
		if (0 == mype && verbosityActivated && bbIterationCounter%1==0) {
			double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);
			BBSMPS_ALG_LOG_SEV(warning)<<"\n----------------------------------------------------\n"<<
			"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":GAP:"<<gap<<":Tree Size:"<<heap.size()<<":Time:"<<BBSMPSSolver::instance()->getWallTime()<<"\n"<<
			"----------------------------------------------------";
		}
	}

//if (0 == mype) BBSMPS_ALG_LOG_SEV(summary) << "Objective function value = " << objUB ;
//if (0 == mype) BBSMPS_ALG_LOG_SEV(summary) << "Objective function LB = " << objLB ;
	if (0 == mype && verbosityActivated) {
		double gap = fabs(objUB-objLB)*100/(fabs(objUB)+10e-10);

		BBSMPS_ALG_LOG_SEV(warning)<<"\n--------------EXPLORATION TERMINATED----------------\n"<<
		"Iteration "<<bbIterationCounter<<":LB:"<<objLB<<":UB:"<<objUB<<":GAP:"<<gap<<":Tree Size:"<<heap.size()<<"\n"<<
		":Nodes Fathomed:"<<nodesFathomed<<":Nodes with integer Solution:"<<nodesBecameInteger<<"\n"<<
		"LP Relaxation Value:"<<LPRelaxationValue<<":LP Relaxation Time:"<<LPRelaxationTime<<":Preprocessing Time:"<<PreProcessingTime<<":Total Time:"<<BBSMPSSolver::instance()->getWallTime();

		"----------------------------------------------------";
		heuristicsManager.printStatistics();
		branchingRuleManager.printStatistics();
		cuttingPlanesManager.printStatistics();
		BBSMPSSolver::instance()->printSolutionStatistics(objLB);
		BBSMPSSolver::instance()->printPresolveStatistics( );
	}
	double t = BBSMPSSolver::instance()->getWallTime();
	if (0 == mype && verbosityActivated) {
		BBSMPS_APP_LOG_SEV(warning)<<boost::format("Branch and Bound took %f seconds") % t;
	}
}


void BBSMPSTree::setTimeLimit(int _tiLim){
	tiLim=_tiLim;
}

void BBSMPSTree::setNodeLimit(int _nodeLim){
	nodeLim=_nodeLim;
}



  void BBSMPSTree::loadSimpleHeuristics(){


	/*BBSMPSHeuristicRounding *hr= new BBSMPSHeuristicRounding(0,1,"SimpleRounding");
   heuristicsManager.addLPHeuristic(hr);

    BBSMPSHeuristicFixAndDive *hr2= new BBSMPSHeuristicFixAndDive(0,1,"FixAndDive");
	heuristicsManager.addLPHeuristic(hr2);

  	BBSMPSHeuristicLockRounding *hr3= new BBSMPSHeuristicLockRounding(0,1,"LockRounding");
  	heuristicsManager.addLPHeuristic(hr3);


*/

  	BBSMPSHeuristicLockRounding *hr6= new BBSMPSHeuristicLockRounding(0,3,"HeuristicLockRounding");

  heuristicsManager.addLPHeuristic(hr6);
/*
    BBSMPSHeuristicFixAndDiveLocks *hr5= new BBSMPSHeuristicFixAndDiveLocks(0,1,"FixAndDiveLocks");

   heuristicsManager.addLPHeuristic(hr5);  */
  }

  void BBSMPSTree::loadCuttingPlanes(){
	BBSMPSCuttingPlaneGenerator01KP *plane=new BBSMPSCuttingPlaneGenerator01KP("01kp");
    cuttingPlanesManager.addCuttingPlaneGenerator(plane);
  }
  void BBSMPSTree::loadMIPHeuristics(){
   	BBSMPSHeuristicRINS *hr= new BBSMPSHeuristicRINS(30,50,"RINS",200);
	BBSMPSHeuristicRENS *hr2= new BBSMPSHeuristicRENS(30,50,"RENS",200);
  /*	BBSMPSHeuristicCrossover *hr3= new BBSMPSHeuristicCrossover(0,1,"Crossover",1);
  	heuristicsManager.addMIPHeuristic(hr3);*/
  heuristicsManager.addMIPHeuristic(hr);
   heuristicsManager.addMIPHeuristic(hr2);
  /*
     BBSMPSHeuristicSolutionRINS *hr4= new BBSMPSHeuristicSolutionRINS(0,1,"SolRINS",1);
   	heuristicsManager.addMIPHeuristic(hr4);


	 BBSMPSHeuristicBestRINSJump *hr5= new BBSMPSHeuristicBestRINSJump(0,1,"BestRinsJump",1);
  	heuristicsManager.addMIPHeuristic(hr5);
    BBSMPSHeuristicSolutionPolishing *hr6= new BBSMPSHeuristicSolutionPolishing(0,1,"BBSMPSHeuristicSolutionPolishing",1);
  	heuristicsManager.addMIPHeuristic(hr6);*/

  }


  void BBSMPSTree::loadLPHeuristic(BBSMPSHeuristic *heur){
  	heuristicsManager.addLPHeuristic(heur);
  }
   void BBSMPSTree::loadMIPHeuristic(BBSMPSHeuristic *heur){
  	heuristicsManager.addMIPHeuristic(heur);
  }


  BBSMPSNode* BBSMPSTree::topOfHeap(){
  	return heap.top();
  }


  void BBSMPSTree::setVerbosity(bool verbose){
  	verbosityActivated=verbose;
  }


void BBSMPSTree::setSolLimit(int _solLim){
	solsDiscoveredInit=BBSMPSSolver::instance()->getSolPoolSize();

   solsDiscoveredLimit=_solLim;

}

void BBSMPSTree::setGAPTolLimit( double _GAPTolLim){
	optGapTol=_GAPTolLim;
}
