#include "BBSMPSHeuristicsManager.hpp"
using namespace std;

BBSMPSHeuristicsManager::BBSMPSHeuristicsManager(){};
BBSMPSHeuristicsManager::~BBSMPSHeuristicsManager(){};


void BBSMPSHeuristicsManager::addHeuristic(BBSMPSHeuristic *heuristic){
	assert(heuristic!=NULL);
	heuristicsList.push_back(heuristic);
}

bool BBSMPSHeuristicsManager::runHeuristics(BBSMPSNode* node,denseBAVector &LPRelaxationSolution, vector<BBSMPSSolution> &solutions, double objUB){
	bool success=false;
	int mype=BBSMPSSolver::instance()->getMype();
	
	for (int it=0;it<heuristicsList.size(); it++){
		BBSMPSHeuristic *heur=heuristicsList[it];
		if (heur->checkPeriodicity(node) && heur->shouldItRun(node,LPRelaxationSolution)){
			BBSMPSSolution auxSol;
			
			bool heurSuccess=heur->runHeuristic(node, LPRelaxationSolution,  auxSol,objUB);
			if (heurSuccess){

				solutions.push_back(auxSol);

			}
			success= success || heurSuccess;
			
		}
		
	}
	return success;
}

void BBSMPSHeuristicsManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(info)<<"++++++++++++++HEURISTIC STATISTICS++++++++++++++++";
	for (int it=0;it<heuristicsList.size(); it++){
		BBSMPSHeuristic *heur=heuristicsList[it];
		heur->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(info)<<"++++++++++++++++++++++++++++++++++++++++++++++++++";
}

void BBSMPSHeuristicsManager::freeResources(){
	for (int it=0;it<heuristicsList.size(); it++){
		BBSMPSHeuristic *heur=heuristicsList[it];
		delete heur;
	}
}