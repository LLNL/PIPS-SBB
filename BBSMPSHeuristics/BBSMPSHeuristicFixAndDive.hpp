// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicRounding.hpp

   Description: Simple nearest integer rounding heuristic

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
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution,double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:

};


#endif
