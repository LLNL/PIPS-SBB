// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicCrossover.hpp

   Description: Solution Crossover PRIMAL HEURISTIC: The solutions in the solution pool 
   are compared two by two. In case we find two solutions that are highly similar, a new 
   subproblem is generated, where the common variables are fixed.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICCROSSOVER_H
#define BBSMPSHEURISTICCROSSOVER_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"
#include <utility> 
#include <map>

class BBSMPSHeuristicCrossover: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicCrossover(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim),BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
	std::map< std::pair<int,int> ,int> seenCrossovers;
};


#endif

