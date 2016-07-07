// ----------------------------------------------------------------------------
/**
   File: BBSMPSHeuristicRENS.hpp

   Description: RENS PRIMAL HEURISTIC: The integral part of the solution that has integer values
   is fixed. Meanwhile, the bounds of the rest of the integer variables are tightened to the nearest
   integer.

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSHEURISTICRENS_H
#define BBSMPSHEURISTICRENS_H

#include "BBSMPSHeuristic.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "BBSMPSTree.hpp"

class BBSMPSHeuristicRENS: public BBSMPSHeuristic {

public:
	BBSMPSHeuristicRENS(int offset, int depth,  const char *_name, int _nodeLim): nodeLim(_nodeLim),BBSMPSHeuristic(offset,depth,_name){};
	bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	int nodeLim;
};


#endif

