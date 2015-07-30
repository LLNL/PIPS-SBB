// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicMagic.hpp

   Description: Buggy version of Timo Berthold's rounding that seemed to do pretty well.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICMAGIC_H
#define BBSMPSHEURISTICMAGIC_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <utility>
class BBSMPSHeuristicMagic: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicMagic(int offset, int depth,  const char *_name): BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution,double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:

};


#endif
