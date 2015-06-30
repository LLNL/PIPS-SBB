#include "BBSMPSSolution.hpp"

	BBSMPSSolution::BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue,double _timeOfDiscovery):
	objValue(_objValue),timeOfDiscovery(_timeOfDiscovery){
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    	solutionVector.allocate(dimsSlacks, ctx, PrimalVector);
		setSolutionVector(_solutionVector);

	}

	BBSMPSSolution::BBSMPSSolution(){
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    	solutionVector.allocate(dimsSlacks, ctx, PrimalVector);
		timeOfDiscovery=-1;

	}

	BBSMPSSolution::~BBSMPSSolution(){}

	void BBSMPSSolution::setSolutionVector(const denseBAVector &_solutionVector){
		solutionVector.copyFrom(_solutionVector);
	}
	void BBSMPSSolution::getSolutionVector(denseBAVector &_solutionVector){
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    	_solutionVector.allocate(dimsSlacks, ctx, PrimalVector);
		_solutionVector.copyFrom(solutionVector);
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

