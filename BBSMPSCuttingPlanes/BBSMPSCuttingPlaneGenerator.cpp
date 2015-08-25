#include "BBSMPSCuttingPlaneGenerator.hpp"

using namespace std;

BBSMPSCuttingPlaneGenerator::BBSMPSCuttingPlaneGenerator(const char *_name){
	timesCalled=0;
	name=_name;
	timesSuccessful=0;

}
BBSMPSCuttingPlaneGenerator::~BBSMPSCuttingPlaneGenerator(){};


void BBSMPSCuttingPlaneGenerator::printStatistics(){
	BBSMPS_ALG_LOG_SEV(info)<<"Heuristic:"<<name<<":Times Called:"<<timesCalled<<":Times successful:"<<timesSuccessful;
}