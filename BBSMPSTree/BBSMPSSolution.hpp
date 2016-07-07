// ----------------------------------------------------------------------------
/**
   File: BBSMPSSolution.hpp

   Description: Class containing all the information relative to a single solution.
   				It encapsulates information such as the actual solution vector as well
   				as the objective value. It will be useful to expand the class to include
   				metadata such as number of fractional variables.

   Limitations: --

*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSSOLUTION_H
#define BBSMPSSOLUTION_H

#include "BAData.hpp"


class BBSMPSSolution {
public:
	bool operator==(const BBSMPSSolution &other) const;
	BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue, double _timeOfDiscovery=-1);
	BBSMPSSolution(){};
	~BBSMPSSolution();
	void setSolutionVector(const denseBAVector &_solutionVector);
	void getSolutionVector(denseBAVector &_solutionVector)const;
	void setObjValue(double _objValue);
	double getObjValue() const;
	double getTimeOfDiscovery() const;
	void setTimeOfDiscovery(const double _timeOfDiscovery);
	int getSolNumber()const;
private:
	//Class variable used to assign solution numbers upon instantiation
	static int solCounter;
	int solNumber;
	denseBAVector solutionVector;
	double objValue;
	double timeOfDiscovery;

};

#endif
