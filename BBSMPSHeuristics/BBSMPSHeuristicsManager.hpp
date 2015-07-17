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
#include "BBSMPSSolver.hpp"

class BBSMPSHeuristicsManager {

public:
	BBSMPSHeuristicsManager();

	~BBSMPSHeuristicsManager();

	void addHeuristic(BBSMPSHeuristic *heuristic);

	bool runHeuristics(BBSMPSNode* n,denseBAVector &LPRelaxationSolution, std::vector<BBSMPSSolution> &solutions, double objUB);

	void printStatistics();

	void freeResources();

private:
	std::vector<BBSMPSHeuristic*> heuristicsList;
};

#endif