// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicsManager.hpp

   Description: Heuristics manager responsible for organizing heuristic rules
   				and calling them recursively upon call.

   Limitations: Heuristics are only called when their frequency allows it or ShouldItRun()
   				returns positive.

*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSHEURISTICSMANAGER_H
#define BBSMPSHEURISTICSMANAGER_H

#include <set>
#include <vector>
#include "BBSMPSHeuristic.hpp"


class BBSMPSHeuristicsManager {

public:
	BBSMPSHeuristicsManager();

	~BBSMPSHeuristicsManager();

	void addHeuristic(BBSMPSHeuristic *heuristic);

	bool runHeuristics(BBSMPSNode* n,denseBAVector &LPRelaxationSolution, std::vector<denseBAVector> &solutions, std::vector<double> &solutionObjVals);

	void printStatistics();

private:
	std::multiset<BBSMPSHeuristic*> heuristicsList;
};

#endif