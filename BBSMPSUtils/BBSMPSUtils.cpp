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
