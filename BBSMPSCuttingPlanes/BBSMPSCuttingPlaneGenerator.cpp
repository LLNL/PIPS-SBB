#include "BBSMPSCuttingPlaneGenerator.hpp"

using namespace std;

BBSMPSCuttingPlaneGenerator::BBSMPSCuttingPlaneGenerator(const char *_name){
	cuttingPlanesGenerated=0;
	name=_name;


}
BBSMPSCuttingPlaneGenerator::~BBSMPSCuttingPlaneGenerator(){};


void BBSMPSCuttingPlaneGenerator::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"Cutting Plane:"<<name<<":Planes Generated:"<<cuttingPlanesGenerated;
}