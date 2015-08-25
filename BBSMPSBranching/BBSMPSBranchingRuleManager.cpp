#include "BBSMPSBranchingRuleManager.hpp"

BBSMPSBranchingRuleManager::BBSMPSBranchingRuleManager(){};
BBSMPSBranchingRuleManager::~BBSMPSBranchingRuleManager(){};


void BBSMPSBranchingRuleManager::addBranchingRule(BBSMPSBranchingRule *rule){
	assert(rule!=NULL);
	branchingRuleList.insert(rule);
}

bool BBSMPSBranchingRuleManager::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln){
	bool success=false;
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end() && !success; ++it){
		BBSMPSBranchingRule *br=(*it);
		success= success || br->branch(node, childNodes,  primalSoln);
		
	}
	return success;
}


void BBSMPSBranchingRuleManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(info)<<"**************HEURISTIC STATISTICS****************";
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end(); ++it){
		BBSMPSBranchingRule *br=(*it);
		br->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(info)<<"**************************************************";
}
