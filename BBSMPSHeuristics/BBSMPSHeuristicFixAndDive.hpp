// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicFixAndDive.hpp

   Description: This diving heuristic proceeds to round one variable at at time, to then proceed to solve the resulting LP relaxation.
   The variable order of choice is the largest fractionality first.
*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICFIXANDDIVE_H
#define BBSMPSHEURISTICFIXANDDIVE_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <utility>
class BBSMPSHeuristicFixAndDive: public BBSMPSHeuristic {

public:
	BBSMPSHeuristicFixAndDive(int offset, int depth,  const char *_name): BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:

};


#endif
