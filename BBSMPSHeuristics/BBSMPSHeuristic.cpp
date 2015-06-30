#include "BBSMPSHeuristic.hpp"

using namespace std;

BBSMPSHeuristic::BBSMPSHeuristic(int _offset, int _depth, const char *_name){
	offset=_offset;
	depth=_depth;
	timesCalled=0;
	name=_name;
	timesSuccessful=0;
	cumulativeTime=0;
}
BBSMPSHeuristic::~BBSMPSHeuristic(){};

bool BBSMPSHeuristic::checkPeriodicity(BBSMPSNode* node){

	int nodeDepth=node->getNodeDepth();
	assert(nodeDepth>=0);
	if (offset>nodeDepth) return false;
	nodeDepth-=offset;
	return (nodeDepth%depth==0);
}


void BBSMPSHeuristic::printStatistics(){
	BBSMPS_ALG_LOG_SEV(summary)<<"Heuristic:"<<name<<":Times Called:"<<timesCalled<<":Times successful:"<<timesSuccessful<<":Execution Time:"<<cumulativeTime;
}