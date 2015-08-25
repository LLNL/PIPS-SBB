// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlaneMIR.hpp

   Description: Mixed Integer Rounding cutting plane

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSCUTTINGPLANEGENERATORMIR_H
#define BBSMPSCUTTINGPLANEGENERATORMIR_H

#include "BBSMPSCuttingPlaneGenerator.hpp"
#include "BBSMPSUtils.hpp"

   
class BBSMPSCuttingPlaneGeneratorMIR: public BBSMPSCuttingPlaneGenerator {
	
public:
	BBSMPSCuttingPlaneGeneratorMIR( const char *_name): BBSMPSCuttingPlaneGenerator(_name){};
	bool generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
};


#endif

