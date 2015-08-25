#include "BBSMPSCuttingPlanesManager.hpp"
using namespace std;

BBSMPSCuttingPlanesManager::BBSMPSCuttingPlanesManager(){};
BBSMPSCuttingPlanesManager::~BBSMPSCuttingPlanesManager(){};


void BBSMPSCuttingPlanesManager::addCuttingPlaneGenerator(BBSMPSCuttingPlaneGenerator *cpg){
	assert(cp!=NULL);
	cuttingPlaneGeneratorsList.push_back(cpg);
}

bool BBSMPSCuttingPlanesManager::generateCuttingPlanes(BBSMPSNode* node,denseBAVector &LPRelaxationSolution){
	bool success=false;
	
	for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){
		BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
		if (cpg->shouldItRun(node,LPRelaxationSolution)){
			bool cutSuccess=cpg->generateCuttingPlane(node, LPRelaxationSolution);
			
			success= success || cutSuccess;
			
		}
		
	}
	return success;
}

void BBSMPSCuttingPlanesManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(info)<<"+++++++++++++CUTTING PLANE STATISTICS+++++++++++++";
	for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){
		BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
		cpg->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(info)<<"++++++++++++++++++++++++++++++++++++++++++++++++++";
}

void BBSMPSCuttingPlanesManager::freeResources(){
	for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){
		BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
		delete cpg;
	}
}