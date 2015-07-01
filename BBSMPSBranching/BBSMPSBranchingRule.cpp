#include "BBSMPSBranchingRule.hpp"

BBSMPSBranchingRule::BBSMPSBranchingRule(int _priority){
	priority=_priority;
    timesCalled=0;
	timesSuccessful=0;
}

BBSMPSBranchingRule::~BBSMPSBranchingRule(){

}

int BBSMPSBranchingRule::getPriority() const{
	return priority;
}

void BBSMPSBranchingRule::setPriority(int _priority){
	priority=_priority;
}

void BBSMPSBranchingRule::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"Branching Rule:"<<name<<":Times Called:"<<timesCalled<<":Times successful:"<<timesSuccessful;
}