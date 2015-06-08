// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicRINS.hpp

   Description: RINS PRIMAL HEURISTIC

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICRINS_H
#define BBSMPSHEURISTICRINS_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"

class BBSMPSHeuristicRINS: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicRINS(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim), BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
};


#endif

