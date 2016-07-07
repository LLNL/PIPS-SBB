#include "BBSMPSCuttingPlanesManager.hpp"
using namespace std;

BBSMPSCuttingPlanesManager::BBSMPSCuttingPlanesManager(int _maxRounds){
	maxRounds=_maxRounds;
};
BBSMPSCuttingPlanesManager::~BBSMPSCuttingPlanesManager(){};


void BBSMPSCuttingPlanesManager::addCuttingPlaneGenerator(BBSMPSCuttingPlaneGenerator *cpg){
	assert(cp!=NULL);
	cuttingPlaneGeneratorsList.push_back(cpg);
}

bool BBSMPSCuttingPlanesManager::generateCuttingPlanes(BBSMPSNode* node,denseBAVector &LPRelaxationSolution){
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	int cuttingPlanesAdded=0;
	denseBAVector primalSoln=denseBAVector(LPRelaxationSolution);
	int nIters=0;

	while(nIters<maxRounds){

		bool success=false;
		rootSolver.setLB(BBSMPSSolver::instance()->getOriginalLB());
		rootSolver.setUB(BBSMPSSolver::instance()->getOriginalUB());

		for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){

			BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
			if (cpg->shouldItRun(node,LPRelaxationSolution)){
				bool cutSuccess=cpg->generateCuttingPlane(node, primalSoln);

				success= success || cutSuccess;

			}

		}
		if (!success) return (cuttingPlanesAdded>0);

		std::vector<BBSMPSCuttingPlane*> newPlanes;
		node->getCurrentNodeCuttingPlanes(newPlanes);
		for (int i=cuttingPlanesAdded; i< newPlanes.size(); i++){
			newPlanes[i]->applyCuttingPlane();
		}
		cuttingPlanesAdded=newPlanes.size();
		BBSMPSSolver::instance()->commitNewColsAndRows();


		denseBAVector lbUpdated(BBSMPSSolver::instance()->getOriginalLB());
		denseBAVector ubUpdated(BBSMPSSolver::instance()->getOriginalUB());

		node->getAllBranchingInformation(lbUpdated,ubUpdated);

		rootSolver.setLB(lbUpdated);
		rootSolver.setUB(ubUpdated);


		BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
		node->reconstructWarmStartState(ps);


		rootSolver.setStates(ps);


		rootSolver.commitStates();

		/* Solve LP defined by current node*/
		//if (0 == mype && verbosityActivated) BBSMPS_ALG_LOG_SEV(summary) << "Solving LP subproblem.";
		rootSolver.go();
		primalSoln=denseBAVector(rootSolver.getPrimalSolution());
		nIters++;
		}


	return (cuttingPlanesAdded>0);

}

void BBSMPSCuttingPlanesManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"=============CUTTING PLANE STATISTICS=============";
	for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){
		BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
		cpg->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(warning)<<"==================================================";
}

void BBSMPSCuttingPlanesManager::freeResources(){
	for (int it=0;it<cuttingPlaneGeneratorsList.size(); it++){
		BBSMPSCuttingPlaneGenerator *cpg=cuttingPlaneGeneratorsList[it];
		delete cpg;
	}
}