//TODO: Figure out more minimal list of includes.
#ifndef PIPS_SBB_PRESOLVE_H
#define PIPS_SBB_PRESOLVE_H

#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include "Presolve.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

// For branch-and-bound code:
#include <queue> // priority queue
#include <cassert> // C-style assertions
#include <cmath> // for floor, ceil, abs functions
#include <algorithm> // for min
#include <limits> // for numeric_limits, overflow checks

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

class Presolve {
public:

  Presolve(BAData& _d, SMPSInput &input) : d(_d),
					   ctx(_d.ctx),
					   dims(input, ctx),
					   dimsSlacks(dims),
					   mype(ctx.mype()),
					   status(LoadedFromFile),
					   intTol(1e-6),
					   eqTol(1e-6),
					   numPresolves(0),
					   numBdsChg(0),
					   numRhsChg(0),
					   numCoeffChg(0),
					   numBinFixed(0),
					   numCtsFixed(0) {

    // Prior to presolve, determine if a given variable & scenario is binary,
    // or integer
    isColInteger.allocate(d.dims, d.ctx, PrimalVector);
    isColBinary.allocate(d.dims, d.ctx, PrimalVector);

    for (int col = 0; col < input.nFirstStageVars(); col++) {
      // Use SMPS input file to determine if variable is integer
      isColInteger.getFirstStageVec()[col] = input.isFirstStageColInteger(col);

      // If variable is not integer, it cannot be binary
      if(!isColInteger.getFirstStageVec()[col]) {
	isColBinary.getFirstStageVec()[col] = false;
      }
      // Otherwise, variable is integer; test bounds to determine if binary
      else {
	double colLB = d.l.getFirstStageVec()[col];
	double colUB = d.u.getFirstStageVec()[col];
	isColBinary.getFirstStageVec()[col] = isBinary(colLB, colUB, intTol);
      }
    }

    // TODO: Fix deep nesting, perhaps by iterating over localScenarios
    // instead of using if(ctx.assignedScenario(scen)) idiom.
    for (int scen = 0; scen < input.nScenarios(); scen++) {
      if(ctx.assignedScenario(scen)) {
	for (int col = 0; col < input.nSecondStageVars(scen); col++) {
	  // Use SMPS input file to determine if variable is integer
	  isColInteger.getSecondStageVec(scen)[col] =
	    input.isSecondStageColInteger(scen, col);

	  // If variable is not integer, it cannot be binary
	  if(!isColInteger.getSecondStageVec(scen)[col]) {
	    isColBinary.getSecondStageVec(scen)[col] = false;
	  }
	  // Otherwise, variable is integer; test bounds to determine if binary
	  else {
	    double colLB = d.l.getSecondStageVec(scen)[col];
	    double colUB = d.u.getSecondStageVec(scen)[col];
	    isColBinary.getSecondStageVec(scen)[col] =
	      isBinary(colLB, colUB, intTol);
	  }
	}
      }
    }

    // Presolve first stage until nothing changes.
    while(true) {
      numPresolves++;
      bool isMIPchanged = false;
      if(0 == mype) PIPS_ALG_LOG_SEV(info) << "First stage presolve iteration "
			 << numPresolves << endl;
      isMIPchanged = presolveFirstStage() || isMIPchanged;

      // Stop presolve if infeasiblility detected after 1st stage presolve.
      if (ProvenInfeasible == status) break;

      if(0 == mype) PIPS_ALG_LOG_SEV(info) << "Second stage presolve iteration "
			 << numPresolves << endl;
      for (int scen = 0; scen < input.nScenarios(); scen++) {
	if(ctx.assignedScenario(scen)) {
	  isMIPchanged = presolveSecondStage(scen) ||
	    isMIPchanged;
	}
      }

      // Synchronize first stage information and isMIPchanged via reductions.
      if(0 == mype) PIPS_ALG_LOG_SEV(info) << "Presolve sync iteration "
			 << numPresolves << endl;
      presolveSyncFirstStage(isMIPchanged);

      // Stop presolve if infeasibility detected after 2nd stage presolve.
      if (ProvenInfeasible == status) break;

      // Stop presolve if presolve operations do not change MIP.
      if (!isMIPchanged) break;

      // Failsafe for now
      if(numPresolves > 100) break;
    }

    logStats();

  }

private:
  // disallow default, copy constructors
  Presolve();
  Presolve(const Presolve& p);
  // disallow copy assignment operator
  Presolve& operator=(const Presolve& p);

  BAContext& ctx;
  BAData &d;

private:
  BADimensions dims;
  BADimensionsSlacks dimsSlacks;
  int mype;

  // TODO: Must centralize these solver settings somewhere
  // Integrality tolerance
  double intTol;

  // Floating point equality tolerance
  double eqTol;

public:
  solverState status;

  // Solver summary statistics
  unsigned int numPresolves, numBdsChg, numRhsChg, numCoeffChg, numBinFixed, numCtsFixed;

private:
  // TODO: Move these fields into something like BAMIPData, etc.
  // Vectors that store whether given variable is binary, integer
  BAFlagVector<bool> isColInteger;
  BAFlagVector<bool> isColBinary;

  void logStats() {
    if(0 == mype) {
      PIPS_ALG_LOG_SEV(info) << "Presolve summary statistics:" << endl;
      PIPS_ALG_LOG_SEV(info) << "Number of presolves: " << numPresolves
			     << endl;
      PIPS_ALG_LOG_SEV(info) << "Number of column bounds changes: " << numBdsChg
			     << endl;
      PIPS_ALG_LOG_SEV(info) << "Number of binary variables fixed: "
			     << numBinFixed << endl;
      PIPS_ALG_LOG_SEV(info) << "Number of row bounds changes: "
			     << numRhsChg << endl;
      PIPS_ALG_LOG_SEV(info) << "Number of coefficient changes: "
			     << numCoeffChg << endl;
    }
  }

  // TODO: Dirty hack. Must refactor status into a class that is a state
  // machine.
  void setStatusToProvenInfeasible() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (ProvenInfeasible == status);
    if (isInReachableState) {
      status = ProvenInfeasible;
    }
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

  // TODO: Move utility functions fracPart & isIntFeas to central location
  // TODO: Eliminate duplication of these functions in PIPS-SBB B&B tree
  // Returns "fractional part" of x
  // Should always be nonnegative.
  double fracPart(double x) {
    return min(x - floor(x), ceil(x) - x);
  }

  // Returns true if "x" is integer feasible up to tolerance "tol"
  bool isIntFeas(double x, double tol) {
    // Alternate method of calculation, useful for benchmarking, unit tests
    //    return ( (abs(floor(x) - x) <= tol) || (abs(ceil(x) - x) <= tol) );
    return (fracPart(x) <= tol);
  }

  double overflowSafeAdd(double x, double y) {
    bool isSameSign = ((x < 0.0) == (y < 0.0));
    bool isMagnitudeOverflow =
      (std::abs(y) > std::numeric_limits<double>::max() - std::abs(x));
    if (isSameSign && isMagnitudeOverflow) {
      if (x > 0.0) return COIN_DBL_MAX;
      if (y < 0.0) return -COIN_DBL_MAX;
    }

    return (x + y);
  }

  // TODO: Refactor presolve into its own class.
  // TODO: Undo changes for bound consistency, check for overflow.

  void incrementLmaxLminByRow(const denseVector &colLB,
			      const denseVector &colUB,
			      const CoinShallowPackedVector &currentRow,
			      double &Lmax,
			      double &Lmin,
			      int scen) {
    const double *ptrToElts = currentRow.getElements();
    const int *ptrToIdx = currentRow.getIndices();
    int currentRowSize = currentRow.getNumElements();

    for (int j = 0; j < currentRowSize; j++) {
      int col = ptrToIdx[j]; // column index from sparse vector
      double coeff = ptrToElts[j];

      // TODO(oxberry1@llnl.gov): Make this code logging code.
      if (0 == mype) {
	PIPS_ALG_LOG_SEV(debug) << "LB of scen " << scen << ", col " << col
				<< " = " << colLB[col] << endl;
	PIPS_ALG_LOG_SEV(debug) << "UB of scen " << scen << ", col " << col
				<< " = " << colUB[col] << endl;
	PIPS_ALG_LOG_SEV(debug) << "coeff of scen " << scen << ", col " << col
				<< " = " << coeff << endl;
      }

      // Here, floating point comparison tolerance is not necessary;
      // only the sign matters.

      if(coeff >= 0) {
	Lmax = overflowSafeAdd(Lmax, colUB[col] * coeff);
	Lmin = overflowSafeAdd(Lmin, colLB[col] * coeff);
      }
      else { // coeff < 0
	//
	//   Note: In Savelsbergh, Section 1.3, first two equations,
	//   coeff is assumed positive by convention, and then a
	//   negative sign is applied to terms with coeff based on
	//   the index j being in a positive index set or negative
	//   index set. Since coeff is negative in this branch, we
	//   flip the sign of those terms, so they are now all
	//   additions instead of subtractions.
	//
	Lmax = overflowSafeAdd(Lmax, colLB[col] * coeff);
	Lmin = overflowSafeAdd(Lmin, colUB[col] * coeff);
      }
    }
  }

  // NOTE: To avoid name clashes, replace "col" by "var" where appropriate,
  // since column and primal variable are synonyms.
  // TODO: Figure out way to keep # of args to 7 or less.
  // TODO: Simplify "bounds change" counter.
  bool tightenColBoundsByRow(denseVector &colLB,
			     denseVector &colUB,
			     const denseFlagVector<bool> &isVarInteger,
			     const denseFlagVector<bool> &isVarBinary,
			     const CoinShallowPackedVector &currentRow,
			     double Lmax,
			     double Lmin,
			     double rowLB,
			     double rowUB,
			     int scen) {
    bool isMIPchanged = false;
    const double *ptrToElts = currentRow.getElements();
    const int *ptrToIdx = currentRow.getIndices();
    int currentRowSize = currentRow.getNumElements();

    // Bound improvement/fixing binary variables/improving binary coeffs
    for (int j = 0; j < currentRowSize; j++) {
      const int col = ptrToIdx[j];
      const bool isInteger = isVarInteger[col];
      const bool isBin = isVarBinary[col]; //TODO: Improve name; isBinary taken
      double coeff = ptrToElts[j];
      bool isCoeffPositive = (coeff > 0); // Note: only nonzero coeffs stored
      bool isCoeffNegative = (coeff < 0);
      double &varUB = colUB[col];
      double &varLB = colLB[col];

      // If coefficient is zero, bounds cannot be improved.
      if (!isCoeffPositive && !isCoeffNegative) continue;

      // If variable bounds are fixed and equal, bounds cannot be improved.
      bool isFixed = (fabs(varUB - varLB) < eqTol);
      if (isFixed) continue;

      // Bound improvement setup steps; note: for valid lower bound, signs
      // of coefficient terms are flipped because Savelsbergh assumes all
      // coefficients are positive in his paper.

      if (isBin) { // binary fixing
	// Use abs value of coefficient because Savelsbergh assumes all
	// coefficients positive.

	// Savelsbergh, Section 1.3: Fixing of variables
	bool isBinaryFixableLmin = (overflowSafeAdd(Lmin, abs(coeff)) > rowUB);
	// Translation of Savelsbergh 1.3: for lower bound
	// constraints, flip "min" to "max", row upper bound to row
	// lower bound, and adding abs val of coefficient to
	// subtracting abs val of coefficient.
	bool isBinaryFixableLmax = (overflowSafeAdd(Lmax, -abs(coeff)) < rowLB);

	// If either binary fixing condition is true, then binary variable is fixable.
	bool isBinaryFixable = isBinaryFixableLmin || isBinaryFixableLmax;

	if(isBinaryFixable) isMIPchanged = true;

	// Four cases:
	// Cases 1 & 2: binary variable fixable based on Lmin
	// and row upper bound
	if (isBinaryFixableLmin) {
	  // Case 1: if coefficient is positive, binary variable fixed to zero
	  if (isCoeffPositive) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Fix bin var in scen "
						   << scen << ", col "
						   << col << " to 0 via Lmin!"
						   << endl;
	    varUB = 0.0;
	    numBdsChg++; numBinFixed++;
	  }
	  // Case 2: if coefficient is negative, binary variable fixed to one
	  if (isCoeffNegative) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Fix bin var in scen "
						   << scen << ", col "
						   << col << " to 1 via Lmin!"
						   << endl;
	    varLB = 1.0;
	    numBdsChg++; numBinFixed++;
	  }
	}

	// Cases 3 & 4: binary variable is fixable based on Lmax and
	// row lower bound
	if (isBinaryFixableLmax) {
	  // Case 3: if coefficient is positive, binary variable fixed to one
	  if (isCoeffPositive) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Fix bin var in scen "
						   << scen << ", col "
						   << col << " to 1 via Lmax!"
						   << endl;
	    varLB = 1.0;
	    numBdsChg++; numBinFixed++;
	  }
	  // Case 4: if coefficient is negative, binary variable fixed to zero
	  if (isCoeffNegative) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Fix bin var in scen "
						   << scen << ", col "
						   << col << " to 0 via Lmax!"
						   << endl;
	    varUB = 0.0;
	    numBdsChg++; numBinFixed++;
	  }
	}

      }
      else { // continuous and integer non-binary variable fixing
	// If assertions in this scope are tripped, check for overflow
	// issues.

	// Lmin-derived bounds -- in Savelsberg, Section 1.1.
	// These expressions both come straight from Savelsbergh's summary
	// in Section 1.3, under "Improvement of bounds".
	if (isCoeffPositive) {
	  double LminUB = (rowUB - (Lmin - coeff*varLB))/coeff;
	  if (LminUB < varUB) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Tightening UB on scen "
						   << scen << ", col " << col
						   << " from " << varUB
						   << " to " << LminUB
						   << "via Lmin!\n";
	    assert(LminUB > varLB);
	    varUB = LminUB;
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	}
	else {
	  double LminLB = (rowUB - (Lmin - coeff*varUB))/coeff;
	  if (LminLB > varLB) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Tightening LB on scen "
						   << scen << ", col " << col
						   << " from " << varLB
						   << " to " << LminLB
						   << " via Lmin!\n";
	    assert(LminLB < varUB);
	    varLB = LminLB;
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	}

	// Lmax-derived bounds -- by analogy to Savelsberg, Section 1.1
	// These expressions both come from Savelsbergh's summary
	// in Section 1.3, under "Improvement of bounds"; the modifications
	// are to flip lower bounds to upper bounds and Lmin to Lmax.
	if (isCoeffPositive) {
	  double LmaxLB = (rowLB - (Lmax - coeff*varUB))/coeff;
	  if (LmaxLB > varLB) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Tightening LB on scen "
						   << scen << ", col " << col
						   << " from " << varLB << " to "
						   << LmaxLB << " via Lmax!\n";
	    assert(LmaxLB < varUB);
	    varLB = LmaxLB;
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	}
	else {
	  double LmaxUB = (rowLB - (Lmax - coeff*varLB))/coeff;
	  if (LmaxUB < varUB) {
	    if (0 == mype) PIPS_ALG_LOG_SEV(debug) << "Tightening UB on scen "
						   << scen << ", col " << col
						   << " from " << varUB
						   << " to " << LmaxUB
						   << " via Lmax!\n";
	    assert(LmaxUB > varLB);
	    varUB = LmaxUB;
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	}

	// TODO: Improve variable bounds by rounding for integer-valued variables
	if(isInteger) {
	  if(!isIntFeas(varLB, intTol)) {
	    varLB = ceil(varLB);
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	  if(!isIntFeas(varUB, intTol)) {
	    varUB = floor(varUB);
	    isMIPchanged = true;
	    numBdsChg++;
	  }
	}
      }
    }
    return isMIPchanged;
  }

  bool improveCoeffsByRow(denseVector &colLB,
			  denseVector &colUB,
			  const denseFlagVector<bool> &isVarBinary,
			  const CoinShallowPackedVector &currentRow,
			  int row,
			  double Lmax,
			  double Lmin,
			  double& rowLB,
			  double& rowUB) {

    bool isMIPchanged = false;

    const double *ptrToElts = currentRow.getElements();
    const int *ptrToIdx = currentRow.getIndices();
    int currentRowSize = currentRow.getNumElements();

    for (int j = 0; j < currentRowSize; j++) {
      const int col = ptrToIdx[j];
      const bool isBin = isVarBinary[col];
      double coeff = ptrToElts[j];
      bool isCoeffPositive = (coeff > 0); // Note: only nonzero coeffs stored
      bool isCoeffNegative = (coeff < 0);
      //double &varUB = colUB[col];
      //double &varLB = colLB[col];

      // Derived by analogy to Savelsbergh Sections 1.2, 1.3
      // Note: This code cannot currently work as written, because
      // PIPSInterface.d (e.g., rootSolver.d) is a protected member, and
      // thus cannot be written to at the moment.
      if(isBin) {
	// Maximum and minimum possible values for inequality in this row
	// after discarding all other inequalities
	double rowMax = overflowSafeAdd(Lmax, -abs(coeff));
	double rowMin = overflowSafeAdd(Lmin,  abs(coeff));

	// Maximum coefficient increases/decreases by coefficient
	// improvement subsection of Savelsbergh, Section 1.2.
	// Note: derived additional relationships for double-sided
	// inequalities. The basic idea is to see how much each side
	// of the inequality could possibly be modified after
	// looking at the max & min possible values for inequality.
	// The minimum of the possible changes will be used to compute
	// changes in bounds and coefficients.
	double coeffLBchg = max(rowUB - rowMax, 0.0);
	double coeffUBchg = max(rowMin - rowLB, 0.0);
	double coeffChg = min(coeffLBchg, coeffUBchg);
	bool isCoeffImprovable = (coeffChg > 0.0);

	// Note: In many cases, if the coefficients cannot be improved
	// (i.e., coeffChg == 0.0), the code below will do unnecessary
	// assignments. One later optimization could be to get rid of
	// these assignments, if they require a significant amount of
	// time.

	// With two-sided inequalities, it is *ALWAYS* possible to
	// improve the coefficient and one side of the bounds (upper
	// or lower).  (Compare to the one-sided inequality case,
	// where it is always possible to improve the coefficient,
	// but not necessarily the bound, which can only be improved
	// for binary variables with positive coefficients.)  The
	// idea is to consider two MIPs at once, and then show that
	// in the positive coefficient case, the minimum of the
	// possible changes can be applied to reduce the row upper
	// bound and the coefficient. In the negative coefficient
	// case, the minimum of the possible changes can be applied
	// to increase the row lower bound and the
	// coefficient. Here, the reason that the negative
	// coefficient is an increase and not a decrease is because
	// Savelsbergh forces coefficients to be positive in his
	// derivations (even the "negative" ones; he implicitly
	// takes an absolute value); flipping signs changes the
	// subtractions to additions.
	// TODO: May need to update isCoeffPositive & isCoeffNegative
	// in response to updates of coeff. Possibly worth putting into
	// its own self-updating data structure?

	if (isCoeffImprovable && isCoeffPositive) {
	  rowUB = overflowSafeAdd(rowUB, -coeffChg);
	  double newCoeff = overflowSafeAdd(coeff, -coeffChg);
	  d.Arow->modifyCoefficient(row, col, newCoeff);
	  // PIPS-S uses the idiom:
	  // Acol->reverseOrderedCopyOf(*Arow);
	  // This idiom seems wasteful here; a possibly better one could be:
	  d.Acol->modifyCoefficient(row, col, newCoeff);
	  numRhsChg++; numCoeffChg++;
	}
	else if (isCoeffImprovable && isCoeffNegative) {
	  rowLB = overflowSafeAdd(rowLB, coeffChg);
	  double newCoeff = overflowSafeAdd(coeff, coeffChg);
	  d.Arow->modifyCoefficient(row, col, newCoeff);
	  // PIPS-S uses the idiom:
	  // Acol->reverseOrderedCopyOf(*Arow);
	  // This idiom seems wasteful here; a possibly better one could be:
	  d.Acol->modifyCoefficient(row, col, newCoeff);
	  numRhsChg++; numCoeffChg++;
	}
      }
    }
    return isMIPchanged;
  }

  // First stage presolve
  bool presolveFirstStage() {

    denseBAVector &lb = d.l;
    denseBAVector &ub = d.u;

    bool isMIPchanged = false;
    if (0 == mype) {
      PIPS_ALG_LOG_SEV(debug) << "First stage has:\n"
			      << "\t " << dims.numFirstStageVars()
			      << " logical variables\n"
			      << "\t " << dims.numFirstStageCons()
			      << " slack variables/constraints\n";
    }

    // Begin presolve.
    // For now, focus only on:
    // - 1st stage variables
    // - upper bound inequalities
    // dims.numFirstStageCons() == dimsSlacks.numFirstStageCons()
    for (int row = 0; row < dims.numFirstStageCons(); row++) {
      // Nomenclature taken from: "Preprocessing and Probing
      // Techniques for Mixed Integer Programming Problems",
      // M. W. P. Savelsbergh, ORSA Journal on Computing,
      // Vol. 6, No. 4, Fall 1994, 445-453.

      // Compute L_{max}^{i} and L_{min}^{i} for row i using Section 1.3
      // of Savelsbergh.
      //
      // Ignoring other constraints, but not bounds:
      // L_max = maximum possible value of constraint in current row
      // L_min = minimum possible value of constraint in current row
      double Lmax = 0, Lmin = 0;
      //      assert(!(d.Arow->isColOrdered()));
      CoinShallowPackedVector currentRow = d.Arow->getVector(row);
      incrementLmaxLminByRow(lb.getFirstStageVec(),
			     ub.getFirstStageVec(),
			     currentRow,
			     Lmax,
			     Lmin,
			     -1);

      // Diagnostic code.
      if (0 == row) {
	if (0 == mype) {
	  PIPS_ALG_LOG_SEV(debug) << "For row 0, Lmin = "
				  << Lmin << " and Lmax = " << Lmax << endl;
	}
      }

      // Constraints are stored in the form l <= Ax <= u. Let nvars =
      // # of variables (including slacks!). For first stage
      // constraint row j, the lower bound on the inequality is stored
      // in lb.getFirstStageVec()[nvars+j], where j goes from 0 to (#
      // of constraints minus 1). The corresponding upper bound on the
      // inequality for first stage constraint row j is stored in
      // ub.getFirstStageVec()[nvars+j].

      // Inferred from BAData::BAData(stochasticInput &input, BAContext &ctx)
      // dimsSlacks.numFirstStageVars() ==
      //    (dims.numFirstStageCons() + dims.numFirstStageVars())
      double &rowUB =
	ub.getFirstStageVec()[dims.numFirstStageVars() + row];
      double &rowLB =
	lb.getFirstStageVec()[dims.numFirstStageVars() + row];

      // Diagnostic code
      if (0 == row) {
	if (0 == mype) {
	  PIPS_ALG_LOG_SEV(debug) << "For row 0, rowUB = " << rowUB
				  << " and rowLB = " << rowLB << endl;
	}
      }

      // If row is infeasible, set status to infeasible, then
      // terminate presolve, noting that MIP has changed.
      bool isRowInfeasible = (Lmin > rowUB) || (Lmax < rowLB);
      if (isRowInfeasible) {
	setStatusToProvenInfeasible();
	if (0 == mype) {
	  PIPS_ALG_LOG_SEV(debug) << "Row " << row
				  << "in Stage 1 is infeasible!" << endl;
	}
	isMIPchanged = true;
	return isMIPchanged;
      }

      // TODO: Add row redundancy check. Not currently implemented
      // because it requires row deletion (to be added).
      /*
      bool isRowRedundant = (Lmax <= rowUB) && (Lmin >= rowLB);
      if (isRowRedundant) {
	// Do something about redundancy; i.e., mark this row for deletion.
	// TODO: Uncomment the line below once row deletion implemented in
	// the body of this if statement.
	// isMIPchanged = true;
	if (0 == mype) {
	    PIPS_ALG_LOG_SEV(debug) << "Row " << row << " is redundant!\n";
	}
	break;
      }
      */

      // Improve bounds and fix binary variables in first stage.
      isMIPchanged = tightenColBoundsByRow(lb.getFirstStageVec(),
					   ub.getFirstStageVec(),
					   isColInteger.getFirstStageVec(),
					   isColBinary.getFirstStageVec(),
					   currentRow,
					   Lmax,
					   Lmin,
					   rowLB,
					   rowUB, -1) || isMIPchanged;

      // Improve coeffs of binary variables in first stage.
      // NOTE: Currently does nothing; requires refactoring PIPSInterface.
      // See function body for details.
      isMIPchanged = improveCoeffsByRow(lb.getFirstStageVec(),
					ub.getFirstStageVec(),
					isColBinary.getFirstStageVec(),
					currentRow,
					row,
					Lmax,
					Lmin,
					rowLB,
					rowUB) || isMIPchanged;
    }
    return isMIPchanged;
  }


  // Second stage presolve
  bool presolveSecondStage(int scen) {

    denseBAVector &lb = d.l;
    denseBAVector &ub = d.u;

    bool isMIPchanged = false;

    if (0 == mype) {
      PIPS_ALG_LOG_SEV(debug) << "Second stage scenario 0 has:\n"
			      << "\t " << dims.numSecondStageVars(0)
			      << " logical variables\n"
			      << "\t " << dims.numSecondStageCons(0)
			      << " slack variables/constraints\n";
    }
    // Only makes sense when called by process that owns scenario scen.
    // NOTE: Probably could replace hard failure with returning
    // isMIPchanged = false with the current logic used for presolves,
    // but this situation could change as the architecture of presolves
    // changes.
    assert(ctx.assignedScenario(scen));

    // Begin second stage presolve
    for (int row = 0; row < dims.numSecondStageCons(scen); row++) {
      // Nomenclature taken from: "Preprocessing and Probing
      // Techniques for Mixed Integer Programming Problems",
      // M. W. P. Savelsbergh, ORSA Journal on Computing,
      // Vol. 6, No. 4, Fall 1994, 445-453.

      // Compute L_{max}^{i} and L_{min}^{i} for row i using Section 1.3
      // of Savelsbergh.
      //
      // Ignoring other constraints, but not bounds:
      // L_max = maximum possible value of constraint in current row
      // L_min = minimum possible value of constraint in current row
      double Lmax = 0, Lmin = 0;
      CoinShallowPackedVector currentTrow, currentWrow;
      currentTrow = d.Trow[scen]->getVector(row);
      currentWrow = d.Wrow[scen]->getVector(row);

      // Increment Lmax & Lmin using row of T matrix
      incrementLmaxLminByRow(lb.getFirstStageVec(),
			     ub.getFirstStageVec(),
			     currentTrow,
			     Lmax,
			     Lmin,
			     -1);

      // Increment Lmax & Lmin using row of W matrix
      incrementLmaxLminByRow(lb.getSecondStageVec(scen),
			     ub.getSecondStageVec(scen),
			     currentWrow,
			     Lmax,
			     Lmin,
			     scen);

      // Diagnostic code.
      if (0 == row) {
	if (0 == mype) {
	  PIPS_ALG_LOG_SEV(debug) << "For row 0, Lmin = " << Lmin
				  << " and Lmax = " << Lmax << endl;
	}
      }

      // Constraints are stored in the form l <= Ax <= u. Let nvars =
      // # of variables (including slacks!). For first stage
      // constraint row j, the lower bound on the inequality is stored
      // in lb.getFirstStageVec()[nvars+j], where j goes from 0 to (#
      // of constraints minus 1). The corresponding upper bound on the
      // inequality for first stage constraint row j is stored in
      // ub.getFirstStageVec()[nvars+j].

      // Inferred from BAData::BAData(stochasticInput &input, BAContext &ctx)
      // dimsSlacks.numFirstStageVars() ==
      //    (dims.numFirstStageCons() + dims.numFirstStageVars())
      double rowUB =
	ub.getSecondStageVec(scen)[dims.numSecondStageVars(scen) + row];
      double rowLB =
	lb.getSecondStageVec(scen)[dims.numSecondStageVars(scen) + row];

      // If row is infeasible, set status to infeasible, then
      // terminate presolve, noting that MIP has changed.
      bool isRowInfeasible = (Lmin > rowUB) || (Lmax < rowLB);
      if (isRowInfeasible) {
	setStatusToProvenInfeasible();
	if (0 == mype) {
	  PIPS_ALG_LOG_SEV(debug) << "Row " << row
				  << " in scenario " << scen
				  << " is infeasible!" << endl;
	}
	isMIPchanged = true;
	return isMIPchanged;
      }

      // TODO: Add row redundancy check. Not currently implemented
      // because it requires row deletion (to be added).
      /*
      bool isRowRedundant = (Lmax <= rowUB) && (Lmin >= rowLB);
      if (isRowRedundant) {
	// Do something about redundancy; i.e., mark this row for deletion.
	// TODO: Uncomment the line below once row deletion implemented in
	// the body of this if statement.
	// isMIPchanged = true;
	if (0 == mype) {
	    PIPS_ALG_LOG_SEV(debug) << "Row " << row << " is redundant!\n";
	}
	break;
      }
      */

      // Tighten first stage column bounds using row of T matrix
      isMIPchanged = tightenColBoundsByRow(lb.getFirstStageVec(),
					   ub.getFirstStageVec(),
					   isColInteger.getFirstStageVec(),
					   isColBinary.getFirstStageVec(),
					   currentTrow,
					   Lmax,
					   Lmin,
					   rowLB,
					   rowUB,
					   -1) || isMIPchanged;

      // Tighten second stage column bounds using row of W matrix
      isMIPchanged = tightenColBoundsByRow(lb.getSecondStageVec(scen),
					   ub.getSecondStageVec(scen),
					   isColInteger.getSecondStageVec(scen),
					   isColBinary.getSecondStageVec(scen),
					   currentWrow,
					   Lmax,
					   Lmin,
					   rowLB,
					   rowUB,
					   scen) || isMIPchanged;

    }

    return isMIPchanged;
  }

  void presolveSyncFirstStage(bool &isMIPchanged) {

    int errorFlag = 0;

    // Detect infeasibilities and broadcast.
    if(ProvenInfeasible == status) {
      errorFlag = MPI_Bcast(&status, 1, MPI_INT, mype, ctx.comm());
    }

    // Note: Must separate this if statement from previous one because
    // the first if statement does communication; if one rank detects
    // infeasibility, it must be broadcast to all ranks. Then we test
    // again on all ranks to return. If the return statement is combined
    // into the previous if statement, there will be a bug because we
    // will only return early on ranks that detect infeasibilities
    // prior to broadcast, which is not the behavior we want.
    if(ProvenInfeasible == status) return;

    denseBAVector &lb = d.l;
    denseBAVector &ub = d.u;

    // TODO: Make the synchronization step more efficient by coalescing
    // communication even more. For now, shoehorn in a less performant
    // implementation as a proof-of-concept. This implementation is the
    // fastest we can do with off-the-shelf reductions; a more performant
    // implementation might implement a custom reduction operation along
    // with a custom buffer used for packing data.
    // TODO: Replace MPI_INT with MPI_LOGICAL all over the place.

    // Min-reduce first stage upper bounds over all ranks
    errorFlag = MPI_Allreduce(MPI_IN_PLACE,
			      ub.getFirstStageVec().getPointer(),
			      dimsSlacks.numFirstStageVars(),
			      MPI_DOUBLE,
			      MPI_MIN,
			      ctx.comm());

    // Max-reduce first stage lower bounds over all ranks
    errorFlag = MPI_Allreduce(MPI_IN_PLACE,
				  lb.getFirstStageVec().getPointer(),
				  dimsSlacks.numFirstStageVars(),
				  MPI_DOUBLE,
				  MPI_MAX,
				  ctx.comm());

    // Implicitly convert bool to int
    int isChanged = isMIPchanged;
    // Logical-OR-reduce isMIPchanged over all ranks
    errorFlag = MPI_Allreduce(MPI_IN_PLACE,
				  &isChanged,
				  1,
				  MPI_INT,
				  MPI_LOR,
				  ctx.comm());
    // NOTE: if a compiler is being a pain about warnings, just negate twice.
    isMIPchanged = static_cast<bool>(isChanged);



  }

};

#endif /* PIPS_SBB_PRESOLVE_H */
