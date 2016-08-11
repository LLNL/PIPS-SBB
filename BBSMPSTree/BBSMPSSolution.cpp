/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
#include "BBSMPSSolution.hpp"

using namespace std;
	int BBSMPSSolution::solCounter=0;

	bool BBSMPSSolution::operator==(const BBSMPSSolution &other) const {
	    return (objValue==other.getObjValue());
	  }

	BBSMPSSolution::BBSMPSSolution(const denseBAVector &_solutionVector, double _objValue,double _timeOfDiscovery):
	objValue(_objValue),timeOfDiscovery(_timeOfDiscovery),solutionVector(_solutionVector){
		solNumber=(++solCounter);
	}



	BBSMPSSolution::~BBSMPSSolution(){}

	void BBSMPSSolution::setSolutionVector(const denseBAVector &_solutionVector){
		solutionVector.copyFrom(_solutionVector);
	}
	void BBSMPSSolution::getSolutionVector(denseBAVector &_solutionVector)const{

		_solutionVector = solutionVector;
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
