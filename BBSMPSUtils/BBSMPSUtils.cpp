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
#include "BBSMPSUtils.hpp"


using namespace std;
// Returns true if "x" is integer feasible up to tolerance "tol"
bool isIntFeas(double x, double tol) {
	return ( (abs(floor(x) - x) <= tol) || (abs(ceil(x) - x) <= tol) );
}

double fracPart(double x) {
	return min(x - floor(x), ceil(x) - x);
}

double roundToNearestInteger(double x){
	if (x-floor(x) < 0.5) return floor(x);
	return ceil(x);
}

double floorFracPart(double x) {
	return x - floor(x);
}

int getFirstStageMinIntInfeasCol(const denseBAVector& primalSoln) {

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	int col;

		    // Return first index of integer variable with fractional value
	for (col = 0; col < input.nFirstStageVars(); col++)
	{

		bool isColInteger = input.isFirstStageColInteger(col);
		bool isValInteger = isIntFeas(primalSoln.getFirstStageVec()[col],
			intTol);

			// If the (col)th 1st stage primal variable is integer,
			// but has a fractional value, return idx
		if(isColInteger && !isValInteger) return col;
	}

		    // Otherwise, 1st stage is integer feasible: return -1;
	return -1;
}




		// NOTE: MPI standard requires passing ints, not bools
int isFirstStageIntFeas(const denseBAVector& primalSoln) {
	return (getFirstStageMinIntInfeasCol(primalSoln) == -1);
}

int getSecondStageMinIntInfeasCol(const denseBAVector& primalSoln, int scen) {
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    // Only makes sense when called by process that owns scenario scen
	assert(ctx.assignedScenario(scen));

	int col;
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col],
			intTol);

			// If the (col)th 2nd stage primal variable of the (scen)th
			// scenario is integer, but has fractional value, return idx
		if (isColInteger && !isValInteger) return col;
	}

		    // Otherwise, return -1;
	return -1;
}

int isSecondStageIntFeas(const denseBAVector& primalSoln, int scen) {
	return (getSecondStageMinIntInfeasCol(primalSoln, scen) == -1);
}


bool isLPIntFeas(const denseBAVector& primalSoln) {

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int is1stStageIntFeas(isFirstStageIntFeas(primalSoln));

		    // Check to see if all 2nd stage scenarios on current MPI rank are
		    // integer feasible. At first scenario that is integer infeasible,
		    // we know that the primal solution is not integer feasible.
	int isMyRankIntFeas(1);
	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			if(!isSecondStageIntFeas(primalSoln, scen)) {
				isMyRankIntFeas = 0;
				break;
			}
		}
	}

	int is2ndStageIntFeas(0);
	int errorFlag = MPI_Allreduce(&isMyRankIntFeas,
		&is2ndStageIntFeas,
		1,
		MPI_INT,
						  MPI_LAND, // MPI logical and
						  ctx.comm());
		    // TODO: Some error handling here

	return (is1stStageIntFeas && is2ndStageIntFeas);
}

// TODO: isZero, isOne, and isBinary are utility methods that should be
  // migrated to a utilities class/namespace. (These functions were given
  // 2 args for flexibility and ease of refactoring.)
  bool isZero(double x, double tol) {
    return (fabs(x) <= tol);
  }

  bool isOne(double x, double tol) {
    return isZero(1 - x, tol);
  }

  bool isBinary(double colLB, double colUB, double tol) {
    return (isZero(colLB, tol) && isOne(colUB, tol));
  }


