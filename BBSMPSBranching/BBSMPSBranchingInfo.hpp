// ----------------------------------------------------------------------------
/**
   File: BBSMPSBranchingInfo.hpp

   Description: Container class for storing the information relative to a single branching
   				decision.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSBRANCHINGINFO_H
#define BBSMPSBRANCHINGINFO_H


class BBSMPSBranchingInfo {

public:
	BBSMPSBranchingInfo(int _varNumber, double _bound, char _direction, int _stage);
	BBSMPSBranchingInfo(int _varNumber, double _bound, char _direction, int _stage, int _scenario);
	~BBSMPSBranchingInfo();

	int getVarNumber();
	double getBound();
	char getDirection();
	int getStageNumber();
	int getScenarioNumber();
	
private:
	int varNumber;
	double bound;
	char direction;
	int stage;
	int scenario;
	BBSMPSBranchingInfo();
};

#endif