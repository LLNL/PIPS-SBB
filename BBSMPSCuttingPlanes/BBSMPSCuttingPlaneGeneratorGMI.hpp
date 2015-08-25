// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlaneGeneratorGMI.hpp

   Description: Gomory mixed integer cutting plane generator

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSCUTTINGPLANEGENERATORGMI_H
#define BBSMPSCUTTINGPLANEGENERATORGMI_H

#include "BBSMPSCuttingPlaneGenerator.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BAData.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSCuttingPlaneGeneratorGMI: public BBSMPSCuttingPlaneGenerator {
	
public:
	BBSMPSCuttingPlaneGeneratorGMI(const char *_name): BBSMPSCuttingPlaneGenerator(_name){};
	bool generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
	bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution);
};


#endif

