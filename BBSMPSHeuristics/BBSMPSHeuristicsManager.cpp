#include "BBSMPSHeuristicsManager.hpp"
using namespace std;

BBSMPSHeuristicsManager::BBSMPSHeuristicsManager(){

};
BBSMPSHeuristicsManager::~BBSMPSHeuristicsManager(){

};

//TODO: Add guards to tell if heuristics need MIP
void BBSMPSHeuristicsManager::addLPHeuristic(BBSMPSHeuristic *heuristic){
	assert(heuristic!=NULL);
	LPHeuristicsList.push_back(heuristic);
}


bool BBSMPSHeuristicsManager::runLPHeuristics(BBSMPSNode* node,denseBAVector &LPRelaxationSolution){
	bool success=false;
	int mype=BBSMPSSolver::instance()->getMype();

	for (int it=0;it<LPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=LPHeuristicsList[it];
		if (heur->checkPeriodicity(node) && heur->shouldItRun(node,LPRelaxationSolution)){
			BBSMPSSolution auxSol;

			bool heurSuccess=heur->runHeuristic(node, LPRelaxationSolution);

			success= success || heurSuccess;

		}

	}

	return success;
}

void BBSMPSHeuristicsManager::addMIPHeuristic(BBSMPSHeuristic *heuristic){
	assert(heuristic!=NULL);
	MIPHeuristicsList.push_back(heuristic);
}

bool BBSMPSHeuristicsManager::runMIPHeuristics(BBSMPSNode* node,denseBAVector &LPRelaxationSolution){
	bool success=false;
	bool first=true;
	for (int it=0;it<MIPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=MIPHeuristicsList[it];
		if (heur->checkPeriodicity(node) && heur->shouldItRun(node,LPRelaxationSolution)){
			BBSMPSSolution auxSol;
			if (first) {
				BBSMPSSolver::instance()->resetSolver();
			}
			first=false;
			bool heurSuccess=heur->runHeuristic(node, LPRelaxationSolution);
			success= success || heurSuccess;

		}

	}
	if (!first)BBSMPSSolver::instance()->resetSolver();
	return success;
}


bool BBSMPSHeuristicsManager::willMIPHeuristicsRun(BBSMPSNode* node,denseBAVector &LPRelaxationSolution){

	bool first=true;
	for (int it=0;it<MIPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=MIPHeuristicsList[it];
		if (heur->checkPeriodicity(node) && heur->shouldItRun(node,LPRelaxationSolution)){
			return true;
		}
	}
	return false;
}


void BBSMPSHeuristicsManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"++++++++++++++HEURISTIC STATISTICS++++++++++++++++";
	double totalHeuristicTime=0;
	for (int it=0;it<LPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=LPHeuristicsList[it];
		heur->printStatistics();
		totalHeuristicTime+=heur->getCumulativeTime();
	}
	for (int it=0;it<MIPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=MIPHeuristicsList[it];
		heur->printStatistics();
		totalHeuristicTime+=heur->getCumulativeTime();
	}
	BBSMPS_ALG_LOG_SEV(warning)<<"Total Heuristic Cumulative Time:"<<totalHeuristicTime;
	BBSMPS_ALG_LOG_SEV(warning)<<"++++++++++++++++++++++++++++++++++++++++++++++++++";
}

void BBSMPSHeuristicsManager::freeResources(){
	for (int it=0;it<LPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=LPHeuristicsList[it];
		if(heur)delete heur;
	}
	for (int it=0;it<MIPHeuristicsList.size(); it++){
		BBSMPSHeuristic *heur=MIPHeuristicsList[it];
		if(heur)delete heur;
	}
}