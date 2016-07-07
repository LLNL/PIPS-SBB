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

	void addLPHeuristic(BBSMPSHeuristic *heuristic);

	void addMIPHeuristic(BBSMPSHeuristic *heuristic);

	bool runLPHeuristics(BBSMPSNode* n,denseBAVector &LPRelaxationSolution);

	bool runMIPHeuristics(BBSMPSNode* n,denseBAVector &LPRelaxationSolution);

	bool willMIPHeuristicsRun(BBSMPSNode* n,denseBAVector &LPRelaxationSolution);

	void printStatistics();

	void freeResources();

private:
	std::vector<BBSMPSHeuristic*> MIPHeuristicsList;
	std::vector<BBSMPSHeuristic*> LPHeuristicsList;

};

#endif