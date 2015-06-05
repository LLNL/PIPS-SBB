// ----------------------------------------------------------------------------
/**
   File: BBSMPSSolver.hpp

   Description: Singleton class that manages the solver and the problem model. It
                contains:
                - The PIPSSInterface
                - The SMPSInput related to the problem model
                - The MPI context
                - The original lower and upper bounds for the problem variables
                - Problem dimensions and slacks

   Limitations: When accessing the solver data, please use the references with caution.
                To instantiate and initialize:

                BBSMPSSolver::initialize(smps); (Where smps is a SMPSInput object)

                Accessing examples:
                
                PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
                BAContext &ctx=BBSMPSSolver::instance()->getBAContext();

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSSOLVER_H
#define BBSMPSSOLVER_H

#include <vector>

#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include "BBSMPSNode.hpp"


class BBSMPSSolver {

public:

  BAContext& getBAContext();
  int& getMype();
  SMPSInput& getSMPSInput();
  PIPSSInterface& getPIPSInterface();
  BADimensions& getBADimensions();
  BADimensionsSlacks& getBADimensionsSlacks();
  denseBAVector& getOriginalLB();
  denseBAVector& getOriginalUB();

  static BBSMPSSolver *instance();
  static BBSMPSSolver *initialize(const SMPSInput &_input);
  static bool isInitialized();

protected:

private:
  BBSMPSSolver(const SMPSInput &_input);

  BAContext ctx; // MPI communication context for PIPS-S
  int mype; // MPI rank of process storing tree (relative to comm in ctx)
  SMPSInput input; // SMPS input file for reading in block angular MILP
  PIPSSInterface rootSolver; // PIPS-S instance for root LP relaxation
  BADimensions dims; // Dimension object for instantiating warm start information
  BADimensionsSlacks dimsSlacks; // Dimension object for warm start info
  denseBAVector lb;
  denseBAVector ub;

  static BBSMPSSolver *solverInstance;
};

#endif