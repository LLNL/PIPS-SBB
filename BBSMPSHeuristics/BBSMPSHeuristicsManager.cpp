#include "BBSMPSHeuristicsManager.hpp"
using namespace std;

BBSMPSHeuristicsManager::BBSMPSHeuristicsManager(){};
BBSMPSHeuristicsManager::~BBSMPSHeuristicsManager(){};


void BBSMPSHeuristicsManager::addHeuristic(BBSMPSHeuristic *heuristic){
	assert(heuristic!=NULL);
	heuristicsList.insert(heuristic);
}

bool BBSMPSHeuristicsManager::runHeuristics(BBSMPSNode* node,denseBAVector &LPRelaxationSolution, vector<denseBAVector> &solutions){
	bool success=false;
	std::multiset<BBSMPSHeuristic*>::iterator it;
	for (it=heuristicsList.begin(); it!=heuristicsList.end(); ++it){
		BBSMPSHeuristic *heur=(*it);
		cout<<"checking iteration "<<endl;
		if (heur->checkPeriodicity(node) && heur->shouldItRun(node,LPRelaxationSolution)){
			denseBAVector aux;
			bool heurSuccess=heur->runHeuristic(node, LPRelaxationSolution,  aux);
			if (heurSuccess){
				solutions.push_back(aux);
			}
			success= success || heurSuccess;
			
		}
		
	}
	return success;
}

void BBSMPSHeuristicsManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(info)<<"++++++++++++++HEURISTIC STATISTICS++++++++++++++++";
	std::multiset<BBSMPSHeuristic*>::iterator it;
	for (it=heuristicsList.begin(); it!=heuristicsList.end(); ++it){
		BBSMPSHeuristic *heur=(*it);
		heur->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(info)<<"++++++++++++++++++++++++++++++++++++++++++++++++++";
}
