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
#include "BBSMPSSolver.hpp"

class BBSMPSSolution {
public:
	BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue, double _timeOfDiscovery=-1);
	BBSMPSSolution();
	~BBSMPSSolution();
	void setSolutionVector(const denseBAVector &_solutionVector);
	void getSolutionVector(denseBAVector &_solutionVector);
	void setObjValue(double _objValue);
	double getObjValue() const;
	double getTimeOfDiscovery() const;
	void setTimeOfDiscovery(const double _timeOfDiscovery);
private:
	denseBAVector solutionVector;
	double objValue;
	double timeOfDiscovery;

};

#endif
