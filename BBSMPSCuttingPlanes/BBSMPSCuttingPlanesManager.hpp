// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlanesManager.hpp

   Description: Cutting planes manager responsible for organizing all the cutting planes
   				and calling them recursively upon application.

   Limitations: Cutting planes are only applied when their frequency allows it and ShouldItRun()
   				returns positive.

*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSCUTTINGPLANESMANAGER_H
#define BBSMPSCUTTINGPLANESMANAGER_H

#include <set>
#include <vector>
#include "BBSMPSCuttingPlaneGenerator.hpp"


class BBSMPSCuttingPlanesManager {

public:
	BBSMPSCuttingPlanesManager();

	~BBSMPSCuttingPlanesManager();

	void addCuttingPlaneGenerator(BBSMPSCuttingPlaneGenerator *cuttingPlaneGen);

	bool generateCuttingPlanes(BBSMPSNode* n,denseBAVector &LPRelaxationSolution);

	void printStatistics();

	void freeResources();

private:
	std::vector<BBSMPSCuttingPlaneGenerator*> cuttingPlaneGeneratorsList;
};

#endif