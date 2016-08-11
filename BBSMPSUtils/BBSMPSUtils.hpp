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
   File: BBSMPSUtils.hpp

   Description: Compendium of auxiliar functions that may be useful/used across
   				multiple classes.

   Limitations: Functions are stateless.

*/
// ----------------------------------------------------------------------------


#ifndef BBSMPSUTILS_H
#define BBSMPSUTILS_H
#include <cmath> // for floor, ceil, abs functions
#include <algorithm> // for min
#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "BBSMPSSolver.hpp"

// Returns true if "x" is integer feasible up to tolerance "tol"
extern bool isIntFeas(double x, double tol);

extern double fracPart(double x);

extern double roundToNearestInteger(double x);

extern double intTol;

  // NOTE: MPI standard requires passing ints, not bools
extern  int isFirstStageIntFeas(const denseBAVector& primalSoln);

extern  int getSecondStageMinIntInfeasCol(const denseBAVector& primalSoln, int scen);

extern  int isSecondStageIntFeas(const denseBAVector& primalSoln, int scen) ;

extern  bool isLPIntFeas(const denseBAVector& primalSoln);

extern double floorFracPart(double x);

// TODO: isZero, isOne, and isBinary are utility methods that should be
  // migrated to a utilities class/namespace. (These functions were given
  // 2 args for flexibility and ease of refactoring.)
extern  bool isZero(double x, double tol);
extern  bool isOne(double x, double tol);
extern  bool isBinary(double colLB, double colUB, double tol);


#endif
