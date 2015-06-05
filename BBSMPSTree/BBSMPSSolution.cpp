#include "BBSMPSSolution.hpp"

	BBSMPSSolution::BBSMPSSolution(denseBAVector &_solutionVector, double _objValue):
	objValue(_objValue){
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    	solutionVector.allocate(dimsSlacks, ctx, PrimalVector);
		setSolutionVector(_solutionVector);

	}

	BBSMPSSolution::BBSMPSSolution(){
		BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    	solutionVector.allocate(dimsSlacks, ctx, PrimalVector);
		

	}

	BBSMPSSolution::~BBSMPSSolution(){}

	void BBSMPSSolution::setSolutionVector(denseBAVector &_solutionVector){
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
	double BBSMPSSolution::getObjValue(){
		return objValue;
	}