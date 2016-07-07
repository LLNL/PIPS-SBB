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
#include "BBSMPSUtils.hpp"

class BBSMPSCuttingPlane {
public:
	BBSMPSCuttingPlane(double lb, double ub, denseBAVector &expr);
	BBSMPSCuttingPlane(double lb, double ub, sparseBAVector &expr);
	BBSMPSCuttingPlane(const BBSMPSCuttingPlane& p);
	BBSMPSCuttingPlane();
	~BBSMPSCuttingPlane();
	bool applyCuttingPlane();
	int getPlaneUid(){return uid;}

private:
	bool applyCrossScenarioCuttingPlane();
	bool applySingleScenarioCuttingPlane();

protected:

		  //Class variable used to assign node numbers upon instantiation
  		static int planeCounter;

		double lb;
		double ub;
		denseBAVector dExpr;
		sparseBAVector sExpr;
		bool isExpressionDense;
		int uid;

};


#endif

