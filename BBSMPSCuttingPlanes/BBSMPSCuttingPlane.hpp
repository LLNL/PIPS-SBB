// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlane.hpp

   Description: Base virtual class of a cutting plane.

   Limitations: An object of this class should not be instantiated, it should always
   be subclassed. 

*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSCUTTINGPLANE_H
#define BBSMPSCUTTINGPLANE_H

#include "BAData.hpp"
#include <cassert> // C-style assertions
#include <string>
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolver.hpp"


class BBSMPSCuttingPlane {
public:
	BBSMPSCuttingPlane(double lb, double ub, denseBAVector &expr);
	~BBSMPSCuttingPlane();
	bool applyCuttingPlane();

private:
	
		
protected:
		
		double lb;
		double ub;
		denseBAVector expr;

};


#endif

