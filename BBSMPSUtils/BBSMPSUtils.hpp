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


#endif