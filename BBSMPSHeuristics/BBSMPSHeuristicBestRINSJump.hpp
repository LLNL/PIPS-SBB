// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicBestRINSJump.hpp

   Description: Solution Crossover PRIMAL HEURISTIC that uses RINS: The best solution in the solution pool 
   is selected and compared against the relaxation. In case both solutions are highly similar, a new 
   subproblem is generated, where the common variables are fixed. We proceed to solve the associated subproblem
   until a better solution is found or a node limit is hit. If a better solution is found, set is as the best solution
   and start the procedure over again.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICBESTRINSJUMP_H
#define BBSMPSHEURISTICBESTRINSJUMP_H

#include "BBSMPSHeuristicLockRounding.hpp"
#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"
#include <utility> 
#include <map>

class BBSMPSHeuristicBestRINSJump: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicBestRINSJump(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim),BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
	std::map< int ,int> seenCrossovers;
};


#endif

