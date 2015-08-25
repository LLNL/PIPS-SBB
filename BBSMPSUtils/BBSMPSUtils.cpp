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
