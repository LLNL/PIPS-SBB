#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include "BBSMPSTree.hpp"
#include "BBSMPSNode.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

// NOTE: Very preliminary code; want to stand up an example
//       first before refining the design

// Use to define branching heuristics?
// enum BranchingRule { FIRST, LAST } ;

// Use to define status of solver?
// Note: "Optimal" conflicts with PIPS-S solverState
//enum BranchAndBoundStatus { Unbounded, Bounded, Feasible, Optimal, Infeasible};

int main(int argc, char **argv) {

        // Initialize MPI
  MPI_Init(&argc, &argv);

        // Get MPI process rank
  int mype;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);

        // Help information if not enough arguments
  if (argc < 2) {
    if (0 == mype) printf("Usage: %s [SMPS root name]\n",argv[0]);
    return 1;
  }

  // Set PIPS logging level.
  
  BBSMPSLogging::init_logging(3);
  //PIPSLogging::init_logging(1 );

  
        // Get SMPS file name and open SMPS file
  string smpsrootname(argv[1]);

  if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Reading SMPS input.";
  SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

  //scoped_ptr<SMPSInput> s(new SMPSInput(datarootname,nscen));

        // Pass communicator to block angular data structures for data distribution
  BAContext ctx(MPI_COMM_WORLD);

  // Initialize branch-and-bound tree
  if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Initializing branch-and-bound tree.";
  BBSMPSTree bb(input);

  bb.loadSimpleHeuristics();
  bb.loadMIPHeuristics();
  bb.setTimeLimit(3600*15);
  //bb.setNodeLimit(5000);
  if (0 == mype) BBSMPS_ALG_LOG_SEV(info) <<"Calling branch-and-bound.";
  bb.branchAndBound();


  if (0 == mype) BBSMPS_APP_LOG_SEV(info) <<"Application successfully terminated.";
  // Clean up MPI data structures
  MPI_Finalize();

  return 0;
}

