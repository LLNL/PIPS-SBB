// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicSolutionPolishing.hpp

   Description: Given the best integer solution, variables are ordered by contribution to
   the objective. The 90% least contributing ones, are fixed. The resulting subproblem is solved
   up to a certain node limit.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICSOLUTIONPOLISHING_H
#define BBSMPSHEURISTICSOLUTIONPOLISHING_H

#include "BBSMPSHeuristicLockRounding.hpp"
#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"
#include <utility> 
#include <map>
#include <vector>
#include <algorithm>
class BBSMPSHeuristicSolutionPolishing: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicSolutionPolishing(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim),BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
	std::map< int ,int> seenCrossovers;
};


#endif

