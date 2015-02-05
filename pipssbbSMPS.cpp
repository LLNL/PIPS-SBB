#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

// For branch-and-bound code:
#include <queue> // priority queue
#include <cassert> // C-style assertions
#include <cmath> // for floor, ceil, abs functions

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

// NOTE: Very preliminary code; want to stand up an example
//       first before refining the design

// Use to define branching heuristics?
// enum BranchingRule { FIRST, LAST } ;

// Use to define status of solver?
// Note: "Optimal" conflicts with PIPS-S solverState
//enum BranchAndBoundStatus { Unbounded, Bounded, Feasible, Optimal, Infeasible};

// Returns true if "x" is integer feasible up to tolerance "tol"
bool isIntFeas(double x, double tol) {
    return ( (abs(floor(x) - x) <= tol) || (abs(ceil(x) - x) <= tol) );
}

class BranchAndBoundNode {
public:
  //PIPSSInterface* solver; // PIPS-S instance
  double parentObj; // obj fn val of node parent

  // TODO: Figure out correct data type for these variables
  denseBAVector lb; // lower bounds of variables
  denseBAVector ub; // upper bounds of variables

  // variable states for warm start information; each index is
  // one of {Basic, AtLower, AtUpper}
  BAFlagVector<variableState> parentStates;

  // TODO: Local cut objects; Global cuts are stored in the B&B tree.

  // Construct node from {lower, upper} bounds, obj val of parent node
  BranchAndBoundNode(double parentObjLB,
		     const denseBAVector &lowerBounds,
		     const denseBAVector &upperBounds,
		     const BAFlagVector<variableState> &states):
    parentObj(parentObjLB),
    lb(lowerBounds),
    ub(upperBounds),
    parentStates(states)  {}

  // Add copy constructor so that priority_queue can use it for
  // instantiation because denseBAVector and BAFlagVector do not have
  // assignment operators or copy constructors.
  BranchAndBoundNode(const BranchAndBoundNode &sourceNode):
    parentObj(sourceNode.parentObj),
    lb(sourceNode.lb),
    ub(sourceNode.ub),
    parentStates(sourceNode.parentStates) {}

  // Overload assignment operator so that priority_queue can use it
  // for instantiation because denseBAVector and BAFlagVector do not
  // have assignment operators or copy constructors. Satisfies
  // "Rule of 3". (For C++11, follow "Rule of 5".)
  BranchAndBoundNode& operator=(const BranchAndBoundNode& sourceNode) {
    // Check for self-assignment
    if (this == &sourceNode) return *this;

    // Copy members:
    parentObj = sourceNode.parentObj;
    lb.copyFrom(sourceNode.lb);
    ub.copyFrom(sourceNode.ub);
    parentStates.copyFrom(sourceNode.parentStates);

    // Return existing object for chaining.
    return *this;
  }

private:
  BranchAndBoundNode(); // Disallow default constructor

};

// Overload the "less than" operator so that priority_queue can use it
// for heapifying comparisons. Make BranchAndBoundNode a templated
// class when adding additional branching heuristics?
bool operator< (const BranchAndBoundNode& left,
		  const BranchAndBoundNode& right)
{
  return left.parentObj < right.parentObj;
}

// TODO: Check if PIPS-S always minimizes; otherwise, must change logic.

class BranchAndBoundTree {
public:
  BAContext ctx; // MPI communication context for PIPS-S
  int mype; // MPI rank of process storing tree (relative to comm in ctx)
  SMPSInput input; // SMPS input file for reading in block angular MILP
  PIPSSInterface rootSolver; // PIPS-S instance for root LP relaxation
  BADimensions dims; // Dimension object for instantiating

  double objUB; // best upper bound on objective function value
  denseBAVector ubPrimalSolution; // primal solution for best UB on obj
  double objLB; // best lower bound on objective function value

  double intTol; // tolerance on integrality checks
  double optGapTol; // tolerance on optimality gap between UB and LB.

  // max-heap data structure with nodes
  // TODO: Refactor to vector<BranchAndBound> & replace w/ make_heap, push_heap, pop_heap
  std::priority_queue<BranchAndBoundNode> heap;

  // Solver status; can only be in the set {LoadedFromFile, Initialized,
  // PrimalFeasible, Optimal, ProvenUnbounded, ProvenInfeasible, Stopped}
  // because there is no duality theory, and currently, the only interface
  // to the solver has to load problem data from a file as one of the first
  // steps.
  solverState status;

public:
  // Upon constructing the B&B tree:
  // - instantiate MPI communicator context for block angular objects
  // - Get rank of process relative to ctx.comm()
  // - instantiate input object
  // - instantiate PIPS-S for root LP relaxation
  // - instantiate dimensions object for holding problem
  //      dimension information for allocating vectors...
  // - instantiate objective function upper bound to +infty (or a value close to that)
  // - set solver status to "LoadedFromFile" because this interface forces MILP to
  //   be loaded from an SMPS file
  BranchAndBoundTree(const SMPSInput& smps): ctx(MPI_COMM_WORLD),
					     mype(ctx.mype()),
					     input(smps),
					     rootSolver(input, ctx, PIPSSInterface::useDual),
					     dims(input, ctx),
					     objUB(COIN_DBL_MAX), objLB(-COIN_DBL_MAX),
					     intTol(1e-6),
					     optGapTol(1e-6),
					     status(LoadedFromFile)
  {

    if (0 == mype) cout << "Calling B&B tree constructor!\n";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start

    BAFlagVector<variableState> states(dims, ctx, PrimalVector);
    // rootSolver.getStates(states);
    rootSolver.setPrimalTolerance(1e-6);
    rootSolver.setDualTolerance(1e-6);

    // TODO: Debug segfault.
    // Bug hypothesis so far: Since I have not yet called "go" on rootSolver
    // and thus I have not yet solved the LP relaxation, the PIPSSInterface
    // object has not yet allocated memory for lower and upper bounds,
    // and also states. Thus, calling the getters for ANY of these members
    // yields a segfault. This hypothesis has been tested empirically.
    if (0 == mype) cout << "This line is still called!\n";
    heap.push(BranchAndBoundNode(objLB,
				 rootSolver.getLB(),
				 rootSolver.getUB(),
				 states));
    if (0 == mype) cout << "First line not called!\n";

  }

  // Default destructor
  ~BranchAndBoundTree() {}

  // Auxiliary functions for branching
  int getFirstStageMinIntInfeasCol(const denseBAVector& primalSoln) {
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
  int isFirstStageIntFeas(const denseBAVector& primalSoln) {
    return (getFirstStageMinIntInfeasCol(primalSoln) == -1);
  }

  int getSecondStageMinIntInfeasCol(const denseBAVector& primalSoln, int scen) {
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

  int isSecondStageIntFeas(const denseBAVector& primalSoln, int scen) {
    return (getSecondStageMinIntInfeasCol(primalSoln, scen) == -1);
  }

  bool isLPIntFeas(const denseBAVector& primalSoln) {

    int is1stStageIntFeas(isFirstStageIntFeas(primalSoln));

    int isMyScenIntFeas(0);
    for (int scen = 0; scen < input.nScenarios(); scen++) {
      if(ctx.assignedScenario(scen)) {
	isMyScenIntFeas = isSecondStageIntFeas(primalSoln, scen);
      }
    }

    int is2ndStageIntFeas(0);
    int errorFlag = MPI_Allreduce(&isMyScenIntFeas,
				  &is2ndStageIntFeas,
				  1,
				  MPI_INT,
				  MPI_LAND, // MPI logical and
				  ctx.comm());
    // TODO: Some error handling here

    return (is1stStageIntFeas && is2ndStageIntFeas);
  }

  void branchOnFirstStage(const denseBAVector& primalSoln) {

    /* Branching setup */
    double lpObj = rootSolver.getObjective();

    /* Get warm start information from solver for B&B node..*/
    BAFlagVector<variableState> states(dims, ctx, PrimalVector);
    rootSolver.getStates(states);

    // Lower (lb) and upper (ub) bounds for two child nodes:
    // the "floor" child and the "ceiling" child
    denseBAVector lbFloor(rootSolver.getLB()), lbCeil(rootSolver.getLB());
    denseBAVector ubFloor(rootSolver.getUB()), ubCeil(rootSolver.getUB());

    // For now, get minimal index of an integer infeasible variable
    int branchCol = getFirstStageMinIntInfeasCol(primalSoln);
    ubFloor.getFirstStageVec()[branchCol] =
      floor(primalSoln.getFirstStageVec()[branchCol]);
    lbCeil.getFirstStageVec()[branchCol] =
      ceil(primalSoln.getFirstStageVec()[branchCol]);

    BranchAndBoundNode childFloor(lpObj, lbFloor, ubFloor, states);
    BranchAndBoundNode childCeil(lpObj, lbCeil, ubCeil, states);
    heap.push(childFloor);
    heap.push(childCeil);

  }

  void branchOnSecondStage(const denseBAVector& primalSoln) {

    /* Branching setup */
    double lpObj = rootSolver.getObjective();

    // Warm start information
    BAFlagVector<variableState> states(dims, ctx, PrimalVector);
    rootSolver.getStates(states);

    // Lower (lb) and upper (ub) bounds for two child nodes:
    // the "floor" child and the "ceiling" child
    denseBAVector lbFloor(rootSolver.getLB()), lbCeil(rootSolver.getLB());
    denseBAVector ubFloor(rootSolver.getUB()), ubCeil(rootSolver.getUB());

    // For now, find the minimum scenario number such that one of its
    // decision variables is integer infeasible. In that scenario number,
    // get the minimal index of an integer infeasible variable.

    int branchableScen(input.nScenarios() + 1);
    int scen = 0;
    for (; scen < input.nScenarios(); scen++)
      {
	if(ctx.assignedScenario(scen)) {
	  if(!isSecondStageIntFeas(primalSoln, scen)) {
	    branchableScen = scen;
	  }
	}
	}
    int branchScen;
    int errorFlag = MPI_Allreduce(&branchableScen,
				  &branchScen,
				  1,
				  MPI_INT,
				  MPI_MIN,
				  ctx.comm());

    // Then, for that scenario number, get the minimal index of
    // an integer infeasible decision variable, and branch on that column
    if(ctx.assignedScenario(branchScen)) {
      int branchCol = getSecondStageMinIntInfeasCol(primalSoln, scen);
      ubFloor.getSecondStageVec(branchScen)[branchCol] =
	floor(primalSoln.getSecondStageVec(branchScen)[branchCol]);
      lbCeil.getSecondStageVec(branchScen)[branchCol] =
	ceil(primalSoln.getSecondStageVec(branchScen)[branchCol]);

      BranchAndBoundNode childFloor(lpObj, lbFloor, ubFloor, states);
      BranchAndBoundNode childCeil(lpObj, lbCeil, ubCeil, states);
      heap.push(childFloor);
      heap.push(childCeil);
    }
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

  void setStatusToPrimalFeasible() {
    bool isInReachableState = (LoadedFromFile == status);
    //      || (status == Bounded);
    if (isInReachableState) {
      status = PrimalFeasible;
    }
  }

  void setStatusToOptimal() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (PrimalFeasible == status); // || (Bounded == status);
    if (isInReachableState) {
      status = Optimal;
    }
  }

  void setStatusToProvenInfeasible() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (ProvenInfeasible == status);
    if (isInReachableState) {
      status = ProvenInfeasible;
    }
  }

  void setStatusToStopped() {
    bool isInReachableState = (LoadedFromFile == status);
    if (isInReachableState) {
      status = Stopped;
    }
  }

  void branchAndBound() {

    unsigned int nodeNumber = 0;
    if (0 == mype) cout << "Starting branch-and-bound!\n";

    /* While heap not empty and there are still nodes in tree */
    // TODO: Add tolerance on optimality gap, time limit option.
    while (true) {

      if (0 == mype) cout << "Heap is empty!\n";

      /* If heap is empty, update status to Stopped (if possible) and break. */
      if (heap.empty()) {
	setStatusToStopped();

	/* If solver status is primal feasible, and the heap is empty, then
	   the solution must be optimal. */
	if (PrimalFeasible == status) {
	  setStatusToOptimal();
	  if (0 == mype) cout << "Optimal solution found!\n";
	}

	/* If solver status is not primal feasible, then the MILP must be
	   infeasible or unbounded. At the moment, there are no checks for
	   boundedness or unboundedness, so return infeasible. */
	// TODO: Add test for unboundedness.
	if (LoadedFromFile == status) {
	  setStatusToProvenInfeasible();
	  if (0 == mype) cout << "MILP is infeasible!\n";
	}

	break;
      }

      /* Get top-most node and pop it off of heap. */
      BranchAndBoundNode currentNode(heap.top());
      heap.pop();
      if (0 == mype) cout << "Popping node " << ++nodeNumber << " off tree!\n";

      /* Set bounds of LP decision variables from BranchAndBoundNode */
      if (0 == mype) cout << "Setting bounds for LP subproblem!\n";
      rootSolver.setLB(currentNode.lb);
      rootSolver.setUB(currentNode.ub);

      /* Set information on basic/nonbasic variables for warm starting */
      if (0 == mype) cout << "Setting warm start information!\n";
      rootSolver.setStates(currentNode.parentStates);
      rootSolver.commitStates();

      /* Solve LP defined by current node*/
      if (0 == mype) cout << "Solving LP subproblem!\n";
      rootSolver.go();

      /* Check solver status for infeasibility/optimality */
      solverState lpStatus = rootSolver.getStatus();

      // Only realistic solver states upon completion:
      // ProvenInfeasible, Optimal, ProvenUnbounded
      // Other solver states are intermediate states that should not
      // hold upon return from rootSolver.
      bool isLPinfeasible = (ProvenInfeasible == lpStatus);
      bool isLPunbounded = (ProvenUnbounded == lpStatus);
      bool isLPoptimal = (Optimal == lpStatus);
      bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
      assert (isLPother); // Error if not infeasible/unbounded/optimal

      /* Fathom by infeasibility */
      // If LP solver returns infeasibility, fathom node, go to start of loop
      if (isLPinfeasible) {
	if (0 == mype) cout << "Fathoming node "
			    << nodeNumber << " by infeasibility!\n";
	continue;
      }

      // Otherwise, LP is optimal or unbounded.
      // If LP solver returns optimal, then the objective is bounded below.
      // TODO: Change solver status to "Bounded".

      // TODO: Combine the integrality and branching steps later

      /* Fathom by value dominance */
      // LP objective function value is lower bound on the objective
      // function value of the LP derived from any node in the subtree
      // of the B&B tree rooted at the current node, so no feasible
      // solution in that subtree can have a lesser objective function
      // value than the current upper bound on the optimal value of
      // the MILP objective function.
      double lpObj = rootSolver.getObjective();
      if (lpObj >= objUB) {
	if (0 == mype) cout << "Fathoming node "
			    << nodeNumber << " by value dominance!\n";
	continue;
      }

      /* Get primal solution */
      denseBAVector primalSoln(rootSolver.getPrimalSolution());

      /* If primal solution is integral: */
      //  - Update solver status to PrimalFeasible
      //  - Check if upper bound improved
      //  - If so, update current primal solution for that upper bound
      //  - Fathom node, go to start of loop

      if(isLPIntFeas(primalSoln)) {
	if (0 == mype) cout << "Node " << nodeNumber << " is integer feasible!\n";
	setStatusToPrimalFeasible();

	// TODO: Maintain solution pool of best k solutions for k = some small value

	/* Update upper bound if it's less than current best upper bound, and
	   the LP solution is optimal (not unbounded). */
	double newUB = rootSolver.getObjective();
	bool isNewUBbetter = (newUB < objUB);
	if (isLPoptimal && isNewUBbetter) {
	  if (0 == mype) cout << "Updating best upper bound to " << newUB << endl;
	  objUB = rootSolver.getObjective();
	  ubPrimalSolution.copyFrom(primalSoln);
	}

	continue;
      }

      // TODO: Fathom by value dominance in breadth-first fashion?

      /* Otherwise, primal solution is not integral: */
      // - Check to see if lower bound can be updated
      // - Branch

      /* Lower bound update */
      // If the LP solver returns an optimal solution AND that solution is
      // greater than the current best lower bound, update the best lower bound
      // on the objective function value.
      if (isLPoptimal) {
	if (lpObj >= objLB) {
	  if (0 == mype) cout << "Updating best lower bound to " << lpObj << endl;
	  lpObj = objLB;
	}
      }

      /* Optimality gap termination criteria: if gap between best
	 upper bound and best lower bound on objective function is
	 less than the optimality gap, stop the solver and return the
	 feasible solution corresponding to the best upper bound. */
      if (abs(objUB - objLB) <= optGapTol) {
	assert(objLB <= objUB);
	setStatusToOptimal();
	if (0 == mype) cout << "Optimality gap reached! Terminating!\n";
	break;
      }

      /* Branching */
      // Decide which stage to branch on:
      // If first stage decision variables not integer feasible,
      // branch on a first stage variable, go to start of loop
      if(!isFirstStageIntFeas(primalSoln)) {
	if (0 == mype) cout << "Branching on first stage!\n";
	branchOnFirstStage(primalSoln);
	continue;
      }

      // If we get to this point, we know that the first stage decision variables
      // are integer feasible, but the LP solution is not integer feasible, so
      // one of the second stage scenarios must not be integer feasible, and
      // one of the variables in one of those scenarios should be branched on.
      if (0 == mype) cout << "Branching on second stage!\n";
      branchOnSecondStage(primalSoln);
      continue;
    }

  }

  private:

  // Make default constructor impossible to call.
  BranchAndBoundTree();

} ;



int main(int argc, char **argv) {

        // Initialize MPI
	MPI_Init(&argc, &argv);

        // Get MPI process rank
	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

        // Help information if not enough arguments
	if (argc < 2) {
		if (mype == 0) printf("Usage: %s [SMPS root name]\n",argv[0]);
		return 1;
	}

        // Get SMPS file name and open SMPS file
	string smpsrootname(argv[1]);

	if (0 == mype) cout << "Reading SMPS input!\n";
	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

	//scoped_ptr<SMPSInput> s(new SMPSInput(datarootname,nscen));

        // Pass communicator to block angular data structures for data distribution
	BAContext ctx(MPI_COMM_WORLD);

	// Initialize branch-and-bound tree
	if (0 == mype) cout << "Initializing branch-and-bound tree!\n";
	BranchAndBoundTree bb(input);

	/*
	// Solve deterministic LP formulation via dual simplex
	PIPSSInterface solver(input, ctx, PIPSSInterface::useDual);

	if (argc == 5) {
		solver.loadStatus(argv[4]);
	}

	//solver.setDumpFrequency(5000,argv[3]);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();

	// Write solution (if given enough input arguments)
	if (argc >= 4 && argv[3][0] != '-') {
		if (mype == 0) printf("Writing solution\n");
		solver.writeStatus(argv[3]);
		if (mype == 0) printf("Finished writing solution\n");
	}
	*/

	if (0 == mype) cout << "Calling branch-and-bound!\n";
	bb.branchAndBound();


	// Clean up MPI data structures
	MPI_Finalize();

	return 0;
}

