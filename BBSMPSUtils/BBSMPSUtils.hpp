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

// Returns true if "x" is integer feasible up to tolerance "tol"
extern bool isIntFeas(double x, double tol); 

extern double fracPart(double x); 

extern double roundToNearestInteger(double x);

extern double intTol;

#endif