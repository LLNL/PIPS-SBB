/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
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
