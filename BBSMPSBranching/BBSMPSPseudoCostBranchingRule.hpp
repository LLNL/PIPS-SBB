// ----------------------------------------------------------------------------
/**
   File: BBSMPSPseudoCostBranchingRule.hpp

   Description: Branching rule that branches on the maximum fractional variable first. 


*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSPSEUDOCOSTBRANCHINGRULE_H
#define BBSMPSPSEUDOCOSTBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSPseudoCostBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);
   	
   	BBSMPSPseudoCostBranchingRule(int priority);


private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMaxFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype) ;

      void initializeVariable(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, double lpRelaxationObjValue, denseBAVector &lb, denseBAVector &ub, int scen, int col);

      bool performRoundOfFirstStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue);

       bool performRoundOfSecondStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue);

       denseBAVector downPseudoCost;
       denseBAVector upPseudoCost;
       denseBAVector downBranchingHistory;
       denseBAVector upBranchingHistory;
       bool everythingFirstStageInitialized;
       bool everythingSecondStageInitialized;
       int reliabilityFactor;

};

#endif

