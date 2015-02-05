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

  double objUB; // upper bound on objective function value
  denseBAVector ubPrimalSolution; // primal solution for upper bound on obj
  double objLB; // lower bound on objective function value

  double intTol; // tolerance on integrality checks

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
					     status(LoadedFromFile)
  {

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start

    BAFlagVector<variableState> states(dims, ctx, PrimalVector);
    rootSolver.getStates(states);

    heap.push(BranchAndBoundNode(objLB,
				 rootSolver.getLB(),
				 rootSolver.getUB(),
				 states));

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

    /* While heap not empty and there are still nodes in tree */
    // TODO: Add tolerance on optimality gap, time limit option.
    while (true) {

      /* If heap is empty, update status to Stopped (if possible) and break. */
      if (heap.empty()) {
	setStatusToStopped();
	break;
      }

      /* Get top-most node and pop it off of heap. */
      BranchAndBoundNode currentNode(heap.top());
      heap.pop();

      /* Set bounds of LP decision variables from BranchAndBoundNode */
      rootSolver.setLB(currentNode.lb);
      rootSolver.setUB(currentNode.ub);

      /* Set information on basic/nonbasic variables for warm starting */
      rootSolver.setStates(currentNode.parentStates);
      rootSolver.commitStates();

      /* Solve LP defined by current node*/
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
      if (isLPinfeasible) continue;

      // Otherwise, LP is optimal or unbounded
      // TODO: Combine the integrality and branching steps later

      /* Fathom by value dominance */
      // LP objective function value is lower bound on the objective
      // function value of the LP derived from any node in the subtree
      // of the B&B tree rooted at the current node, so no feasible
      // solution in that subtree can have a lesser objective function
      // value than the current upper bound on the optimal value of
      // the MILP objective function.
      // TODO: Check if PIPS-S always minimizes; otherwise, must change logic.
      double lpObj = rootSolver.getObjective();

      // TODO: Add tolerance check
      if (lpObj >= objUB) continue;

      /* Get primal solution */
      denseBAVector primalSoln(rootSolver.getPrimalSolution());

      /* Fathom if integer feasible*/
         
      // If primal solution is integral:
      //  - Update upper bound on objective function value
      //  - Update current primal solution for that upper bound
      //  - Fathom node, go to start of loop

      if(isLPIntFeas(primalSoln)) {
	// TODO: Maintain solution pool of best k solutions for k = some small value

	/* Update upper bound if it's less than current best upper bound. */
	double newUB = rootSolver.getObjective();
	if (newUB < objUB) {
	  objUB = rootSolver.getObjective(); 
	  ubPrimalSolution.copyFrom(primalSoln);
	}
	continue;
      }
      
      /* Otherwise, branch. */
      
      // TODO: Fathom by value dominance in breadth-first fashion?

      // Decide which stage to branch on:
      // If first stage decision variables not integer feasible,
      // branch on a first stage variable, go to start of loop
      if(!isFirstStageIntFeas(primalSoln)) {
	branchOnFirstStage(primalSoln);
	continue;
      }
      
      // If we get to this point, we know that the first stage decision variables
      // are integer feasible, but the LP solution is not integer feasible, so
      // one of the second stage scenarios must not be integer feasible, and
      // one of the variables in one of those scenarios should be branched on.
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

	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

	//scoped_ptr<SMPSInput> s(new SMPSInput(datarootname,nscen));
	
        // Pass communicator to block angular data structures for data distribution
	BAContext ctx(MPI_COMM_WORLD);

	// Initialize branch-and-bound tree



	// Solve deterministic LP formulation via dual simplex
	PIPSSInterface solver(input, ctx, PIPSSInterface::useDual);

	if (argc == 5) {
		solver.loadStatus(argv[4]);
	}

	//solver.setDumpFrequency(5000,argv[3]);	
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();

	// Get solution for warm start


	// (Set solution for warm start?)


	// Change bounds...
	
	
	// Write solution (if given enough input arguments)
	if (argc >= 4 && argv[3][0] != '-') {
		if (mype == 0) printf("Writing solution\n");
		solver.writeStatus(argv[3]);
		if (mype == 0) printf("Finished writing solution\n");
	}

	// Clean up MPI data structures
	MPI_Finalize();

	return 0;
}

