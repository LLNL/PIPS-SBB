// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicSolutionRINS.hpp

   Description: Solution Crossover with the LP relaxation

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICSOLUTIONRINS_H
#define BBSMPSHEURISTICSOLUTIONRINS_H

#include "BBSMPSHeuristicLockRounding.hpp"
#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"
#include <utility>
#include <map>

class BBSMPSHeuristicSolutionRINS: public BBSMPSHeuristic {

public:
	BBSMPSHeuristicSolutionRINS(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim),BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
};


#endif

