// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlane01KP.hpp

   Description: 0-1 Knapsack cover cutting plane

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSCUTTINGPLANEGENERATOR01KP_H
#define BBSMPSCUTTINGPLANEGENERATOR01KP_H

#include "BBSMPSCuttingPlaneGenerator.hpp"
#include "BBSMPSUtils.hpp"
#include <utility>
#include <algorithm>
#include "BBSMPSCuttingPlane.hpp"
#include "BBSMPSUtils.hpp"
class BBSMPSCuttingPlaneGenerator01KP: public BBSMPSCuttingPlaneGenerator {

public:
	BBSMPSCuttingPlaneGenerator01KP( const char *_name);
	bool generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);

private:
	sparseBAVector knapsackRows;
	int totalKnapsackRows;

};


#endif

