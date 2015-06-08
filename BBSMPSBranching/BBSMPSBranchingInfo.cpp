#include "BBSMPSBranchingInfo.hpp"

using namespace std;


BBSMPSBranchingInfo::BBSMPSBranchingInfo(int _varNumber, double _bound, char _direction, int _stage, int _scenario):
varNumber(_varNumber),
bound(_bound),
direction(_direction),
stage(_stage),
scenario(_scenario){}
BBSMPSBranchingInfo::~BBSMPSBranchingInfo(){

}

BBSMPSBranchingInfo::BBSMPSBranchingInfo(){

}

int BBSMPSBranchingInfo::getVarNumber(){
	return varNumber;
}
double BBSMPSBranchingInfo::getBound(){
	return bound;
}
char BBSMPSBranchingInfo::getDirection(){
	return direction;
}

int BBSMPSBranchingInfo::getStageNumber(){
	return stage;
}

int BBSMPSBranchingInfo::getScenarioNumber(){
	return scenario;
}
