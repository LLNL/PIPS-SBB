// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicLockRounding.hpp

   Description: Heuristic lock rounding implemented as described by Timo Berthold.
   

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICLOCKROUNDING_H
#define BBSMPSHEURISTICLOCKROUNDING_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <utility>
#include "CoinShallowPackedVector.hpp"
class BBSMPSHeuristicLockRounding: public BBSMPSHeuristic {
	
public:
	BBSMPSHeuristicLockRounding(int offset, int depth,  const char *_name): BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution,double objUB);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:

};


#endif
