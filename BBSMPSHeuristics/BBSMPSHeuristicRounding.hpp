// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicRounding.hpp

   Description: Simple nearest integer rounding heuristic

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICROUNDING_H
#define BBSMPSHEURISTICROUNDING_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
class BBSMPSHeuristicRounding: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicRounding(int offset, int depth,  const char *_name): BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
};


#endif

