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
#include "Presolve.hpp"
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolution.hpp"
#include <set>

struct solutionComparison {
  bool operator() (const BBSMPSSolution& lhs, const BBSMPSSolution& rhs) const
  {return (lhs.getObjValue()<rhs.getObjValue());}
};

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
  const BAFlagVector<variableState>& getOriginalWarmStart();
  void setOriginalWarmStart(BAFlagVector<variableState>&warmStart);
  denseBAVector& getLPRelaxation();
  void setLPRelaxation(denseBAVector &_LPRelaxation);
  void setLPRelaxationObjectiveValue(double lpRelObjVal);
  double getLPRelaxationObjectiveValue();
  static BBSMPSSolver *instance();
  static BBSMPSSolver *initialize(const SMPSInput &_input);
  static bool isInitialized();
  static void deInitialize();
  void printPresolveStatistics();
  void addSolutionToPool(BBSMPSSolution &sol);
  void printSolutionStatistics(double objLB);
  const BBSMPSSolution &getSoln(int index);
  const BBSMPSSolution &getSolnBySolNumber(int number);
  
  int getSolPoolSize();
  double getWallTime();

protected:

private:
  BBSMPSSolver(const SMPSInput &_input);
  BAContext ctx; // MPI communication context for PIPS-S
  int mype; // MPI rank of process storing tree (relative to comm in ctx)
  SMPSInput input; // SMPS input file for reading in block angular MILP
   BAData problemData;
  PIPSSInterface rootSolver; // PIPS-S instance for root LP relaxation
    Presolve pre;
  BADimensions dims; // Dimension object for instantiating warm start information
  BADimensionsSlacks dimsSlacks; // Dimension object for warm start info
  denseBAVector lb;
  denseBAVector ub;
  denseBAVector LPRelaxation;
  BAFlagVector<variableState> originalWarmStart;
  std::set<BBSMPSSolution,solutionComparison> solutionPool;
  
  double LPRelaxationObjectiveValue;
  double startTimeStamp;
  static BBSMPSSolver *solverInstance;
};

#endif