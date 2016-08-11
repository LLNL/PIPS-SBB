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

