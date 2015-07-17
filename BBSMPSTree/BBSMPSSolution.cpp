#include "BBSMPSSolution.hpp"

using namespace std;
	int BBSMPSSolution::solCounter=0;

	bool BBSMPSSolution::operator==(const BBSMPSSolution &other) const {
	    return (objValue==other.getObjValue());
	  }

	BBSMPSSolution::BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue,double _timeOfDiscovery):
	objValue(_objValue),timeOfDiscovery(_timeOfDiscovery),solutionVector(_solutionVector){
		solNumber=(++solCounter);
		cout<<"Generating solution with "<<solNumber<<" "<<_objValue<<" "<<_timeOfDiscovery<<endl;
	}



	BBSMPSSolution::~BBSMPSSolution(){}

	void BBSMPSSolution::setSolutionVector(const denseBAVector &_solutionVector){
		solutionVector.copyFrom(_solutionVector);
	}
	void BBSMPSSolution::getSolutionVector(denseBAVector &_solutionVector)const{
		
		_solutionVector = denseBAVector(solutionVector);
	}
	void BBSMPSSolution::setObjValue(double _objValue){
		objValue=_objValue;
	}
	double BBSMPSSolution::getObjValue()const{
		return objValue;
	}

	double BBSMPSSolution::getTimeOfDiscovery() const{
		return timeOfDiscovery;
	}

	void BBSMPSSolution::setTimeOfDiscovery(const double _timeOfDiscovery){
	  timeOfDiscovery=_timeOfDiscovery;
	}

	int BBSMPSSolution::getSolNumber()const{
		return solNumber;
	}