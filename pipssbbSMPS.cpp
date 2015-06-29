#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
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

// NOTE: Very preliminary code; want to stand up an example
//       first before refining the design

// Use to define branching heuristics?
// enum BranchingRule { FIRST, LAST } ;

// Use to define status of solver?
// Note: "Optimal" conflicts with PIPS-S solverState
//enum BranchAndBoundStatus { Unbounded, Bounded, Feasible, Optimal, Infeasible};

// Returns true if "x" is integer feasible up to tolerance "tol"
bool isIntFeas(double x, double tol) {
    return ( (abs(floor(x) - x) <= tol) || (abs(ceil(x) - x) <= tol) );
}

double fracPart(double x) {
  return min(x - floor(x), ceil(x) - x);
}

// Outputs solver status:
void outputLPStatus(solverState lpStatus) {
  cout << "PIPS-S has returned ";
  string status;
  switch(lpStatus) {
  case Uninitialized:
    status = "Uninitialized";
    break;
  case LoadedFromFile:
    status = "LoadedFromFile";
    break;
  case Initialized:
    status = "Initialized";
    break;
  case PrimalFeasible:
    status = "PrimalFeasible";
    break;
  case DualFeasible:
    status = "DualFeasible";
    break;
  case Optimal:
    status = "Optimal";
    break;
  case ProvenUnbounded:
    status = "ProvenUnbounded";
    break;
  case ProvenInfeasible:
    status = "ProvenInfeasible";
    break;
  case Stopped:
    status = "Stopped";
    break;
  }
  cout << status << endl;
}

class Presolve {
public:

  Presolve(BAData d, SMPSInput &input) : d(d),
					  input(input),
					  ctx(d.ctx),
					  dims(input, ctx),
					  dimsSlacks(dims),
					  mype(ctx.mype()),
					  status(LoadedFromFile) {

    // Prior to presolve, determine if a given variable & scenario is binary.
    isColBinary.allocate(d.dims, d.ctx, PrimalVector);
    for (int col = 0; col < input.nFirstStageVars(); col++) {
      isColBinary.getFirstStageVec()[col] = input.isFirstStageColInteger(col);
    }
    for (int scen = 0; scen < input.nScenarios(); scen++) {
      if(ctx.assignedScenario(scen)) {
	for (int col = 0; col < input.nSecondStageVars(scen); col++) {
	  isColBinary.getSecondStageVec(scen)[col] =
	    input.isSecondStageColInteger(scen,col);
	}
      }
    }

    // Presolve first stage until nothing changes.
    unsigned int numPresolves = 0;
    while(true) {
      numPresolves++;
      bool isMIPchanged = false;
      if(0 == mype) cout << "First stage presolve iteration "
			 << numPresolves << endl;
      isMIPchanged = presolveFirstStage() || isMIPchanged;

      // Stop presolve if infeasiblility detected after 1st stage presolve.
      if (ProvenInfeasible == status) break;

      if(0 == mype) cout << "Second stage presolve iteration "
			 << numPresolves << endl;
      for (int scen = 0; scen < input.nScenarios(); scen++) {
	if(ctx.assignedScenario(scen)) {
	  isMIPchanged = presolveSecondStage(scen) ||
	    isMIPchanged;
	}
      }

      // Synchronize first stage information and isMIPchanged via reductions.
      if(0 == mype) cout << "Presolve sync iteration "
			 << numPresolves << endl;
      presolveSyncFirstStage(isMIPchanged);

      // Stop presolve if infeasibility detected after 2nd stage presolve.
      if (ProvenInfeasible == status) break;

      // Stop presolve if presolve operations do not change MIP.
      if (!isMIPchanged) break;

      // Failsafe for now
      if(numPresolves > 100) break;
    }

    if (0 == mype) cout << "There were " << numPresolves << " presolves.\n";

  }

private:
  // disallow default, copy constructors
  Presolve();
  Presolve (const Presolve& p); 
  // disallow copy assignment operator
  Presolve& operator=(const Presolve& p);

  SMPSInput input;
  BAContext& ctx;

public:
  BAData d;

private:
  BADimensions dims;
  BADimensionsSlacks dimsSlacks;
  int mype;

public:
  solverState status;

private:
  // TODO: Move this field into something like BAMIPData, etc.
  BAFlagVector<bool> isColBinary; // stores whether a given variable is binary

  // TODO: Dirty hack. Must refactor status into a class that is a state
  // machine.
  void setStatusToProvenInfeasible() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (ProvenInfeasible == status);
    if (isInReachableState) {
      status = ProvenInfeasible;
    }
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
	cout << "LB of scen " << scen << ", col " << col
	     << " = " << colLB[col] << endl;
	cout << "UB of scen " << scen << ", col " << col
	     << " = " << colUB[col] << endl;
	cout << "coeff of scen " << scen << ", col " << col
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
  bool tightenColBoundsByRow(denseVector &colLB,
			     denseVector &colUB,
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
      const bool isBinary = isVarBinary[col];
      double coeff = ptrToElts[j];
      bool isCoeffPositive = (coeff > 0); // Note: only nonzero coeffs stored
      bool isCoeffNegative = (coeff < 0);
      double &varUB = colUB[col];
      double &varLB = colLB[col];

      // If coefficient is zero, bounds cannot be improved.
      if (!isCoeffPositive && !isCoeffNegative) continue;

      // If variable bounds are fixed and equal, bounds cannot be improved.
      bool isFixed = (varUB == varLB);
      if (isFixed) continue;

      // Bound improvement setup steps; note: for valid lower bound, signs
      // of coefficient terms are flipped because Savelsbergh assumes all
      // coefficients are positive in his paper.

      if (isBinary) { // binary fixing
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
	    if (0 == mype) cout << "Fix bin var in scen " << scen << ", col "
				<< col << " to 0 via Lmin!" << endl;
	    varUB = 0.0;
	  }
	  // Case 2: if coefficient is negative, binary variable fixed to one
	  if (isCoeffNegative) {
	    if (0 == mype) cout << "Fix bin var in scen " << scen << ", col "
				<< col << " to 1 via Lmin!" << endl;
	    varLB = 1.0;
	  }
	}

	// Cases 3 & 4: binary variable is fixable based on Lmax and
	// row lower bound
	if (isBinaryFixableLmax) {
	  // Case 3: if coefficient is positive, binary variable fixed to one
	  if (isCoeffPositive) {
	    if (0 == mype) cout << "Fix bin var in scen " << scen << ", col "
				<< col << " to 1 via Lmax!" << endl;
	    varLB = 1.0;
	  }
	  // Case 4: if coefficient is negative, binary variable fixed to zero
	  if (isCoeffNegative) {
	    if (0 == mype) cout << "Fix bin var in scen " << scen << ", col "
				<< col << " to 0 via Lmax!" << endl;
	    varUB = 0.0;
	  }
	}

      }
      else { // continuous and integer non-binary variable fixing
	// If assertions in this scope are tripped, check for overflow
	// issues.

	// Lmin-derived lower bounds -- in Savelsberg, Section 1.1.
	// These expressions both come straight from Savelsbergh's summary
	// in Section 1.3, under "Improvement of bounds".
	if (isCoeffPositive) {
	  double LminUB = (rowUB - (Lmin - coeff*varLB))/coeff;
	  if (LminUB < varUB) {
	    if (0 == mype) cout << "Tightening UB on scen " << scen
				<< ", col " << col << " from "
				<< varUB << " to " << LminUB
				<< "via Lmin!\n";
	    assert(LminUB > varLB);
	    varUB = LminUB;
	    isMIPchanged = true;
	  }
	}
	else {
	  double LminLB = (rowUB - (Lmin - coeff*varUB))/coeff;
	  if (LminLB > varLB) {
	    if (0 == mype) cout << "Tightening LB on scen " << scen
				<< ", col " << col << " from "
				<< varLB << " to " << LminLB
				<< " via Lmin!\n";
	    assert(LminLB < varUB);
	    varLB = LminLB;
	    isMIPchanged = true;
	  }
	}

	// Lmax-derived lower bounds -- by analogy to Savelsberg, Section 1.1
	// These expressions both come from Savelsbergh's summary
	// in Section 1.3, under "Improvement of bounds"; the modifications
	// are to flip lower bounds to upper bounds and Lmin to Lmax.
	if (isCoeffPositive) {
	  double LmaxLB = (rowLB - (Lmax - coeff*varUB))/coeff;
	  if (LmaxLB > varLB) {
	    if (0 == mype) cout << "Tightening LB on scen " << scen
				<< ", col " << col << " from "
				<< varLB << " to " << LmaxLB
				<< " via Lmax!\n";
	    assert(LmaxLB < varUB);
	    varLB = LmaxLB;
	    isMIPchanged = true;
	  }
	}
	else {
	  double LmaxUB = (rowLB - (Lmax - coeff*varLB))/coeff;
	  if (LmaxUB < varUB) {
	    if (0 == mype) cout << "Tightening UB on scen " << scen
				<< ", col " << col << " from "
				<< varUB << " to " << LmaxUB
				<< " via Lmax!\n";
	    assert(LmaxUB > varLB);
	    varUB = LmaxUB;
	    isMIPchanged = true;
	  }
	}

	// TODO: Improve variable bounds by rounding for integer-valued variables
      }
    }
    return isMIPchanged;
  }

  bool improveCoeffsByRow(denseVector &colLB,
			  denseVector &colUB,
			  const denseFlagVector<bool> &isVarBinary,
			  const CoinShallowPackedVector &currentRow,
			  const double &Lmax,
			  const double &Lmin,
			  const double &rowLB,
			  const double &rowUB) {
    bool isMIPchanged = false;

    /*
    const double *ptrToElts = currentRow.getElements();
    const int *ptrToIdx = currentRow.getIndices();
    int currentRowSize = currentRow.getNumElements();

    for (int j = 0; j < currentRowSize; j++) {
      const int col = ptrToIdx[j];
      const bool isBinary = isVarBinary[col];
      double coeff = ptrToElts[j];
      bool isCoeffPositive = (coeff > 0); // Note: only nonzero coeffs stored
      bool isCoeffNegative = (coeff < 0);
      double &varUB = colUB[col];
      double &varLB = colLB[col];

      // Derived by analogy to Savelsbergh Sections 1.2, 1.3
      // Note: This code cannot currently work as written, because
      // PIPSInterface.d (e.g., rootSolver.d) is a protected member, and
      // thus cannot be written to at the moment.
      if(isBinary) {
      // Maximum and minimum possible values for inequality in this row
      // after discarding all other inequalities
      double rowMax = Lmax - abs(coeff);
      double rowMin = Lmin + abs(coeff);

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
      ub.getFirstStageVec()[dimsSlacks.numFirstStageVars() + row] += -coeffChg;
      coeff -= coeffChg;
      problemData.Arow->modifyCoefficient(row, col, newCoeff);
      // PIPS-S uses the idiom:
      // Acol->reverseOrderedCopyOf(*Arow);
      // This idiom seems wasteful here; a possibly better one could be:
      problemData.Acol->modifyCoefficient(row, col, newCoeff);
      }
      else if (isCoeffImprovable && isCoeffNegative) {
      lb.getFirstStageVec()[dimsSlacks.numFirstStageVars() + row] += coeffChg;
      double coeff += coeffChg;
      problemData.Arow->modifyCoefficient(row, col, newCoeff);
      // PIPS-S uses the idiom:
      // Acol->reverseOrderedCopyOf(*Arow);
      // This idiom seems wasteful here; a possibly better one could be:
      problemData.Acol->modifyCoefficient(row, col, newCoeff);
      }
      }
    */
    return isMIPchanged;
  }

  // First stage presolve
  bool presolveFirstStage() {

    denseBAVector &lb = d.l;
    denseBAVector &ub = d.u;

    bool isMIPchanged = false;
    if (0 == mype) {
      cout << "First stage has:\n"
	   << "\t " << dims.numFirstStageVars() << " logical variables\n"
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
	  cout << "For row 0, Lmin = " << Lmin << " and Lmax = " << Lmax << endl;
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
	ub.getFirstStageVec()[dims.numFirstStageVars() + row];
      double rowLB =
	lb.getFirstStageVec()[dims.numFirstStageVars() + row];

      // Diagnostic code
      if (0 == row) {
	if (0 == mype) {
	  cout << "For row 0, rowUB = " << rowUB << " and rowLB = " << rowLB << endl;
	}
      }

      // If row is infeasible, set status to infeasible, then
      // terminate presolve, noting that MIP has changed.
      bool isRowInfeasible = (Lmin > rowUB) || (Lmax < rowLB);
      if (isRowInfeasible) {
	setStatusToProvenInfeasible();
	if (0 == mype) {
	  cout << "Row " << row << "in Stage 1 is infeasible!" << endl;
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
	    cout << "Row " << row << " is redundant!\n";
	}
	break;
      }
      */

      // Improve bounds and fix binary variables in first stage.
      isMIPchanged = tightenColBoundsByRow(lb.getFirstStageVec(),
					   ub.getFirstStageVec(),
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
      cout << "Second stage scenario 0 has:\n"
	   << "\t " << dims.numSecondStageVars(0) << " logical variables\n"
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
	  cout << "For row 0, Lmin = " << Lmin << " and Lmax = " << Lmax << endl;
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
	  cout << "Row " << row
	       << " in scenario " << scen << " is infeasible!" << endl;
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
	    cout << "Row " << row << " is redundant!\n";
	}
	break;
      }
      */

      // Tighten first stage column bounds using row of T matrix
      isMIPchanged = tightenColBoundsByRow(lb.getFirstStageVec(),
					   ub.getFirstStageVec(),
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
    int errorFlag = MPI_Allreduce(MPI_IN_PLACE,
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

    // Detect infeasibilities and broadcast.
    if(ProvenInfeasible == status) {
      errorFlag = MPI_Bcast(&status, 1, MPI_INT, mype, ctx.comm());
    }

  }

};

class BranchAndBoundNode {
public:
  //PIPSSInterface* solver; // PIPS-S instance
  double parentObj; // obj fn val of node parent

  // TODO: Figure out correct data type for these variables
  denseBAVector lb; // lower bounds of variables
  denseBAVector ub; // upper bounds of variables

  // variable states for warm start information; each index is
  // one of {Basic, AtLower, AtUpper}
  BAFlagVector<variableState> parentStates;

  // TODO: Local cut objects; Global cuts are stored in the B&B tree.

  // Construct node from {lower, upper} bounds, obj val of parent node
  BranchAndBoundNode(double parentObjLB,
		     const denseBAVector &lowerBounds,
		     const denseBAVector &upperBounds,
		     const BAFlagVector<variableState> &states):
    parentObj(parentObjLB),
    lb(lowerBounds),
    ub(upperBounds),
    parentStates(states) {}

  // Add copy constructor so that priority_queue can use it for
  // instantiation because denseBAVector and BAFlagVector do not have
  // assignment operators or copy constructors.
  BranchAndBoundNode(const BranchAndBoundNode &sourceNode):
    parentObj(sourceNode.parentObj),
    lb(sourceNode.lb),
    ub(sourceNode.ub),
    parentStates(sourceNode.parentStates) {}

  // Overload assignment operator so that priority_queue can use it
  // for instantiation because denseBAVector and BAFlagVector do not
  // have assignment operators or copy constructors. Satisfies
  // "Rule of 3". (For C++11, follow "Rule of 5".)
  BranchAndBoundNode& operator=(const BranchAndBoundNode& sourceNode) {

    /* Diagnostic code: delete when no longer needed */
    /* Also assumes communicator is MPI_COMM_WORLD. */
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    //    if (0 == mype) cout << "Calling copy assignment operator!\n";
    // Check for self-assignment
    if (this == &sourceNode) {

      //if (0 == mype) cout << "Calling self-assignment branch!\n";
      return *this;
    }

    // Copy-assign each member individually
    parentObj = sourceNode.parentObj;
    lb = sourceNode.lb;
    ub = sourceNode.ub;
    parentStates = sourceNode.parentStates;

    // Return existing object for chaining.
    //if (0 == mype) cout << "Exiting copy assignment operator!\n";
    return *this;
  }

private:
  BranchAndBoundNode(); // Disallow default constructor

};

// Overload the "less than" operator so that priority_queue can use it
// for heapifying comparisons. Make BranchAndBoundNode a templated
// class when adding additional branching heuristics?
bool operator< (const BranchAndBoundNode& left,
		  const BranchAndBoundNode& right)
{
  return left.parentObj < right.parentObj;
}

bool operator> (const BranchAndBoundNode& left,
		  const BranchAndBoundNode& right)
{
  return left.parentObj > right.parentObj;
}


// TODO: Check if PIPS-S always minimizes; otherwise, must change logic.

class BranchAndBoundTree {
public:
  BAContext ctx; // MPI communication context for PIPS-S
  int mype; // MPI rank of process storing tree (relative to comm in ctx)
  SMPSInput input; // SMPS input file for reading in block angular MILP
  BADimensions dims; // Dimension object for instantiating warm start information
  BADimensionsSlacks dimsSlacks; // Dimension object for warm start info
  BAData problemData; // Data structure encoding MIP.
  //  Presolve pre; // Presolve object that transforms LP with presolve ops
  // Solver status; can only be in the set {LoadedFromFile, Initialized,
  // PrimalFeasible, Optimal, ProvenUnbounded, ProvenInfeasible, Stopped}
  // because there is no duality theory, and currently, the only interface
  // to the solver has to load problem data from a file as one of the first
  // steps.
  solverState status;
  PIPSSInterface rootSolver; // PIPS-S instance for root LP relaxation

  double objUB; // best upper bound on objective function value
  denseBAVector ubPrimalSolution; // primal solution for best UB on obj
  double objLB; // best lower bound on objective function value

  double intTol; // tolerance on integrality checks
  double optGapTol; // tolerance on optimality gap between UB and LB.
  double lpPrimalTol; // tolerance on LP primal problems
  double lpDualTol; // tolerance on LP dual problems
  double compTol; // tolerance on LP objective function comparisons

  // max-heap data structure with nodes
  // TODO: Refactor to vector<BranchAndBound> & replace w/ make_heap, push_heap, pop_heap
  //std::priority_queue<BranchAndBoundNode, std::vector<BranchAndBoundNode>, std::less<BranchAndBoundNode> > heap; // max-heap
  std::priority_queue<BranchAndBoundNode, std::vector<BranchAndBoundNode>, std::greater<BranchAndBoundNode> > heap; // min-heap


public:
  // Upon constructing the B&B tree:
  // - instantiate MPI communicator context for block angular objects
  // - Get rank of process relative to ctx.comm()
  // - instantiate input object
  // - instantiate PIPS-S for root LP relaxation
  // - instantiate dimensions object for holding problem
  //      dimension information for allocating vectors...
  // - instantiate objective function upper bound to +infty (or a value close to that)
  // - set integrality tolerance
  // - set optimality gap tolerances
  // - set LP solver tolerances
  // - set comparison tolerance to primal tolerance for now
  // - set solver status to "LoadedFromFile" because this interface forces MILP to
  //   be loaded from an SMPS file

  BranchAndBoundTree(const SMPSInput& smps): ctx(MPI_COMM_WORLD),
					     mype(ctx.mype()),
					     input(smps),
					     problemData(input, ctx),
					     //pre(problemData, input),
					     //status(pre.status),
					     status(LoadedFromFile),
					     //rootSolver(problemData, PIPSSInterface::useDual),
					     rootSolver(input, ctx, PIPSSInterface::useDual),
					     dims(input, ctx),
					     dimsSlacks(dims),
					     objUB(COIN_DBL_MAX),
					     objLB(-COIN_DBL_MAX),
					     intTol(1e-6),
					     optGapTol(1e-6),
					     lpPrimalTol(1e-6),
					     lpDualTol(1e-6),
					     compTol(lpPrimalTol)
  {

    //if (0 == mype) cout << "Calling B&B tree constructor!\n";

    /* Initialize branch-and-bound tree/heap */
    // Get {lower, upper} bounds on decision variables, lower bound on objective function
    // value from parent LP, initialize a node, and push onto heap to start.
    assert (heap.empty()); // heap should be empty to start

    rootSolver.setPrimalTolerance(lpPrimalTol);
    rootSolver.setDualTolerance(lpDualTol);

    if (0 == mype) cout << "Getting problem data from root solver!\n";
    //    BAData d = rootSolver.getBAData();
    if (0 == mype) cout << "BAData object has dimensions:\n"
			<< problemData.dims.numFirstStageVars()
			<< " first stage variables\n"
			<< problemData.dims.numFirstStageCons()
			<< " first stage constraints\n";

    denseBAVector lb, ub;
    if (0 == mype) cout << "Getting lower bounds from root solver!\n";
    lb.allocate(problemData.dims, problemData.ctx, PrimalVector); lb.copyFrom(problemData.l);
    if (0 == mype) cout << "Getting upper bounds from root solver!\n";
    ub.allocate(problemData.dims, problemData.ctx, PrimalVector); ub.copyFrom(problemData.u);
    //    rootSolver.setLB(lb);
    //    rootSolver.setUB(ub);

    // Allocate current best primal solution; normally this primal solution
    // is for the upper bound, but here, we have only the solution to an
    // LP relaxation, which may not be primal feasible. We don't check
    // primal/integer feasibility here.
    //if (0 == mype) cout << "Allocating primal solution!" << endl;
    ubPrimalSolution.allocate(problemData.dims, problemData.ctx, PrimalVector);
    //if (0 == mype) cout << "Getting primal solution!" << endl;
    //ubPrimalSolution.copyFrom(rootSolver.getPrimalSolution());
    //if (0 == mype) cout << "MIP Primal solution updated!" << endl;

    // State of primal variables + slacks for warm starts; slacks must
    // be included because these are used in a reformulation of the
    // problem to standard form
    BAFlagVector<variableState> states(problemData.dims, problemData.ctx, PrimalVector);
    //BAFlagVector<variableState> states;

    // Push root node onto B&B tree/heap.
    BranchAndBoundNode rootNode(objLB, lb, ub, states);
    //if (0 == mype) cout << "Pushing root node onto B&B tree!\n";
    heap.push(rootNode);
    //if (0 == mype) cout << "Exiting B&B constructor!\n";

  }

  // Default destructor
  ~BranchAndBoundTree() {}


  // Auxiliary functions for branching
  int getFirstStageMinIntInfeasCol(const denseBAVector& primalSoln) {
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

  int getFirstStageMaxFracPartCol(const denseBAVector& primalSoln){
    int col;

    double maxFracPart = 0;
    int maxFracPartCol = -1;

    // Return index of integer variable with largest fractional part
    for (col = 0; col < input.nFirstStageVars(); col++)
      {
	bool isColInteger = input.isFirstStageColInteger(col);
	double colFracPart = fracPart(primalSoln.getFirstStageVec()[col]);

	if(isColInteger && (colFracPart > maxFracPart) ) {
	  maxFracPartCol = col;
	  maxFracPart = colFracPart;
	}
      }

    if (maxFracPart <= intTol) return -1;
    return maxFracPartCol;

    }

  // NOTE: MPI standard requires passing ints, not bools
  int isFirstStageIntFeas(const denseBAVector& primalSoln) {
    return (getFirstStageMinIntInfeasCol(primalSoln) == -1);
  }

  int getSecondStageMinIntInfeasCol(const denseBAVector& primalSoln, int scen) {
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

  void branchOnFirstStage(const denseBAVector& primalSoln) {

    /* Branching setup */
    double lpObj = rootSolver.getObjective();

    /* Get warm start information from solver for B&B node..*/
    BAFlagVector<variableState> states(problemData.dims, ctx, PrimalVector);
    rootSolver.getStates(states);

    // Lower (lb) and upper (ub) bounds for two child nodes:
    // the "floor" child and the "ceiling" child
    denseBAVector lbFloor(rootSolver.getLB()), lbCeil(rootSolver.getLB());
    denseBAVector ubFloor(rootSolver.getUB()), ubCeil(rootSolver.getUB());

    /* Branching Rule */
    // For now, get minimal index of an integer infeasible variable
    //int branchCol = getFirstStageMinIntInfeasCol(primalSoln);

    // Get index of maximum fractional part.
    int branchCol = getFirstStageMaxFracPartCol(primalSoln);
    assert(branchCol > -1); // Should always be true if not integer feasible

    if (0 == mype) cout << "Branching on first stage variable "
			<< branchCol << "!\n";
    ubFloor.getFirstStageVec()[branchCol] =
      floor(primalSoln.getFirstStageVec()[branchCol]);
    lbCeil.getFirstStageVec()[branchCol] =
      ceil(primalSoln.getFirstStageVec()[branchCol]);

    BranchAndBoundNode childFloor(lpObj, lbFloor, ubFloor, states);
    BranchAndBoundNode childCeil(lpObj, lbCeil, ubCeil, states);
    heap.push(childFloor);
    heap.push(childCeil);

  }

  void branchOnSecondStage(const denseBAVector& primalSoln) {

    /* Branching setup */
    double lpObj = rootSolver.getObjective();

    // Warm start information
    BAFlagVector<variableState> states(problemData.dims, ctx, PrimalVector);
    rootSolver.getStates(states);

    // Lower (lb) and upper (ub) bounds for two child nodes:
    // the "floor" child and the "ceiling" child
    denseBAVector lbFloor(rootSolver.getLB()), lbCeil(rootSolver.getLB());
    denseBAVector ubFloor(rootSolver.getUB()), ubCeil(rootSolver.getUB());

    // For now, find the minimum scenario number on each rank such
    // that one of its decision variables is integer infeasible.
    // Call that scenario number the branching candidate for each rank.
    // Then find the minimum branching candidate over all ranks. Branch on
    /// that scenario.

    // In the scenario number selected for branching, get the minimal
    // index of an integer infeasible variable.

    int myRankBranchScen(input.nScenarios() + 1);
    //if (0 == mype) cout << "myRankBranchScen = " << myRankBranchScen << endl;
    for (int scen = 0; scen < input.nScenarios(); scen++)
      {
	if(ctx.assignedScenario(scen)) {
	  if(!isSecondStageIntFeas(primalSoln, scen)) {
	    myRankBranchScen = scen;
	    break;
	  }
	}
      }

    int branchScen;
    int errorFlag = MPI_Allreduce(&myRankBranchScen,
				  &branchScen,
				  1,
				  MPI_INT,
				  MPI_MIN,
				  ctx.comm());
    if (0 == mype) cout << "Branching on second stage scenario "
			<< branchScen << "!\n";
    cout << "Processor " << mype << " will branch on second stage scenario "
	 << branchScen << "!\n";

    // Then, for that scenario number, get the minimal index of
    // an integer infeasible decision variable, and branch on that column
    if(ctx.assignedScenario(branchScen)) {
      int branchCol = getSecondStageMinIntInfeasCol(primalSoln, branchScen);
      ubFloor.getSecondStageVec(branchScen)[branchCol] =
	floor(primalSoln.getSecondStageVec(branchScen)[branchCol]);
      lbCeil.getSecondStageVec(branchScen)[branchCol] =
	ceil(primalSoln.getSecondStageVec(branchScen)[branchCol]);
    }

      BranchAndBoundNode childFloor(lpObj, lbFloor, ubFloor, states);
      BranchAndBoundNode childCeil(lpObj, lbCeil, ubCeil, states);
      heap.push(childFloor);
      heap.push(childCeil);
  }

  // TODO: Add way to access single scenario.

  /* TODO: Refactor MIP solver status into its own class to handle
     the transition logic. */
  /* TODO: Add more nuanced version of solver statuses. */
  /* Methods defining allowable transitions for status */
  /* status transition state diagram */
  // NOTE: "Bounded" is not currently a solver state;
  // TODO: Add Bounded as a solver state.
  // LoadedFromFile -> {Bounded, PrimalFeasible, Optimal,
  //                    ProvenInfeasible, Unbounded, Stopped}
  // Bounded -> {PrimalFeasible, Optimal, ProvenInfeasible, Stopped}
  // PrimalFeasible -> {Optimal, Stopped}

  // void setStatusToBounded() {
  //   bool isInReachableState = (LoadedFromFile == status);
  //   if (isInReachableState) {
  //     status = Bounded;
  //   }
  // }

  void setStatusToPrimalFeasible() {
    bool isInReachableState = (LoadedFromFile == status);
    //      || (status == Bounded);
    if (isInReachableState) {
      status = PrimalFeasible;
    }
  }

  void setStatusToOptimal() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (PrimalFeasible == status); // || (Bounded == status);
    if (isInReachableState) {
      status = Optimal;
    }
  }

  void setStatusToProvenInfeasible() {
    bool isInReachableState = (LoadedFromFile == status) ||
      (ProvenInfeasible == status);
    if (isInReachableState) {
      status = ProvenInfeasible;
    }
  }

  void setStatusToStopped() {
    bool isInReachableState = (LoadedFromFile == status);
    if (isInReachableState) {
      status = Stopped;
    }
  }

  void branchAndBound() {

    unsigned int nodeNumber = 0;
    if (0 == mype) cout << "Starting branch-and-bound!\n";

    /* While heap not empty and there are still nodes in tree */
    // TODO: Add tolerance on optimality gap, time limit option.
    while (true) {

      /* If heap is empty, update status to Stopped (if possible) and break. */
      if (heap.empty()) {
	if (0 == mype) cout << "Heap is empty!\n";
	setStatusToStopped();

	/* If solver status is primal feasible, and the heap is empty, then
	   the solution must be optimal. */
	if (PrimalFeasible == status) {
	  setStatusToOptimal();
	  if (0 == mype) cout << "Optimal solution found!\n";
	}

	/* If solver status is not primal feasible, then the MILP must be
	   infeasible. */
	// TODO: Add test for unboundedness.
	if (LoadedFromFile == status) {
	  setStatusToProvenInfeasible();
	  if (0 == mype) cout << "MILP is infeasible!\n";
	}

	break;
      }

      /* Get top-most node and pop it off of heap. */
      if (0 == mype) cout << "Copying node " << ++nodeNumber << " off tree!\n";
      BranchAndBoundNode currentNode(heap.top());
      if (0 == mype) cout << "Popping node " << nodeNumber << " off tree!\n";
      heap.pop();

      if (nodeNumber > 1) { // Otherwise, this node isn't the root node of the B&B tree
	// Set bounds of LP decision variables from BranchAndBoundNode
	if (0 == mype) cout << "Setting bounds for LP subproblem!\n";
	rootSolver.setLB(currentNode.lb);
	rootSolver.setUB(currentNode.ub);
      
	// Set information on basic/nonbasic variables for warm starting
	if (0 == mype) cout << "Setting warm start information!\n";
	// Only set warm start states if not root node;
	// trying to set the warm start states for the root node without
	// a known basic feasible solution will crash PIPS-S.
	rootSolver.setStates(currentNode.parentStates);
	rootSolver.commitStates();
      }

      /* Solve LP defined by current node*/
      if (0 == mype) cout << "Solving LP subproblem!\n";
      rootSolver.go();

      /*
      if (nodeNumber == 1) {
	denseBAVector lb, ub;
	lb.allocate(dims, ctx, PrimalVector);
	ub.allocate(dims, ctx, PrimalVector);
	lb.copyFrom(rootSolver.getBAData().l);
	ub.copyFrom(rootSolver.getBAData().u);
	BAFlagVector<variableState> states(dims, ctx, PrimalVector);
	rootSolver.getStates(states);
	BranchAndBoundNode newCurrentNode(objLB, lb, ub, states);
	currentNode = newCurrentNode;
      }
      */

      /* Check solver status for infeasibility/optimality */
      solverState lpStatus = rootSolver.getStatus();

      // Only realistic solver states upon completion:
      // ProvenInfeasible, Optimal, ProvenUnbounded
      // Other solver states are intermediate states that should not
      // hold upon return from rootSolver.
      if (0 == mype) outputLPStatus(lpStatus);

      bool isLPinfeasible = (ProvenInfeasible == lpStatus);
      //if (0 == mype) cout << "isLPinfeasible = " << isLPinfeasible << endl;
      bool isLPunbounded = (ProvenUnbounded == lpStatus);
      //if (0 == mype) cout << "isLPunbounded = " << isLPunbounded << endl;
      bool isLPoptimal = (Optimal == lpStatus);
      //if (0 == mype) cout << "isLPoptimal = " << isLPoptimal << endl;
      bool isLPother = (!isLPinfeasible && !isLPunbounded && !isLPoptimal);
      //if (0 == mype) cout << "isLPother = " << isLPother << endl;
      assert (!isLPother); // Error if not infeasible/unbounded/optimal

      /* Fathom by infeasibility */
      // If LP solver returns infeasibility, fathom node, go to start of loop
      if (0 == mype) cout << "Checking for infeasibility...\n";
      if (isLPinfeasible) {
	if (0 == mype) cout << "Fathoming node "
			    << nodeNumber << " by infeasibility!\n";
	continue;
      }

      // Otherwise, LP is feasible. LP may be optimal or unbounded.
      // If LP is unbounded, so is the MILP.
      if (isLPunbounded) {
	if (0 == mype) cout << "LP relaxation of node " << nodeNumber
			    << " is unbounded!\n"
			    << "Please add additional constraints to "
			    << "bound the MILP.\n";
	return;
      }

      // At this point, LP must be optimal.
      assert (isLPoptimal); // Error if not optimal.

      // TODO: Combine the integrality and branching steps later

      /* If LP solution is optimal, can fathom by value dominance. */

      if (isLPoptimal) {
      // If LP solver returns optimal, then the objective is bounded below.
      // TODO: Change solver status to "Bounded".
	//	setStatusToBounded();

	/* Lower bound update */ // This part is still incorrect!
	// If the LP solver returns an optimal solution AND that solution is
	// greater than the current best lower bound, update the best lower bound
	// on the objective function value.

	// Since the branch-and-bound tree is stored as a min-heap, the current
	// node being explored always has the minimal objective function value.
	// If the branch-and-bound tree is ever re-heapified so that it is NOT
	// a min-heap, but has the heap property for some other ordering, then
	// this update cannot occur without taking the min objective function
	// value over all values of the parent objective function for nodes
	// still in the B&B tree (which is expensive if it is not the key used
	// to heapify the min-heap).
	//objLB = currentNode.parentObj;
	//if ((lpObj - compTol) >= objLB) {
	//  if (0 == mype) cout << "Current best lower bound is " << objLB << endl;
	//  if (0 == mype) cout << "Updating best lower bound to " << lpObj << endl;
	//  objLB = lpObj;
	//}

	/* Fathom by value dominance */
	// Optimal LP objective function value is lower bound on the objective
	// function value of the LP derived from any node in the subtree
	// of the B&B tree rooted at the current node, so no feasible
	// solution in that subtree can have a lesser objective function
	// value than the current upper bound on the optimal value of
	// the MILP objective function.

	// Get LP objective function value.
	if (0 == mype) cout << "Getting LP objective...\n";
	double lpObj = rootSolver.getObjective();

	if (0 == mype) cout << "Checking for value dominance...\n";
	if ((lpObj - compTol) >= objUB) {
	  if (0 == mype) cout << "Fathoming node "
			      << nodeNumber << " by value dominance!\n";
	  continue;
	}
      }

      /* Get primal solution */
      if (0 == mype) cout << "Getting primal solution...\n";
      denseBAVector primalSoln(rootSolver.getPrimalSolution());

      /* If primal solution is integral: */
      //  - Update solver status to PrimalFeasible
      //  - Check if upper bound improved
      //  - If so, update current primal solution for that upper bound
      //  - Fathom node, go to start of loop

      if (0 == mype) cout << "Checking for integrality of primal solution...\n";
      if(isLPIntFeas(primalSoln)) {
	if (0 == mype) cout << "Node " << nodeNumber << " is integer feasible!\n";
	setStatusToPrimalFeasible();

	// TODO: Maintain solution pool of best k solutions for k = some small value

	/* Update upper bound if it's less than current best upper bound, and
	   the LP solution is optimal (not unbounded). */
	double newUB = rootSolver.getObjective();
	bool isNewUBbetter = (newUB < (objUB - compTol));
	if (isLPoptimal && isNewUBbetter) {
	  if (0 == mype) cout << "Updating best upper bound to " << newUB << endl;
	  objUB = rootSolver.getObjective();
	  ubPrimalSolution.copyFrom(primalSoln);
	}

	continue;
      }

      // TODO: Fathom by value dominance in breadth-first fashion?

      /* Otherwise, primal solution is not integral: */
      // - Check to see if lower bound can be updated
      // - Branch

      /* Optimality gap termination criteria: if gap between best
	 upper bound and best lower bound on objective function is
	 less than the optimality gap, stop the solver and return the
	 feasible solution corresponding to the best upper bound. */
      if (abs(objUB - objLB) <= optGapTol) {
	assert(objLB <= objUB); // this statement could be tripped by numerical error
	setStatusToOptimal();
	if (0 == mype) cout << "Optimality gap reached! Terminating!\n";
	break;
      }

      /* Branching */
      // Decide which stage to branch on:
      // If first stage decision variables not integer feasible,
      // branch on a first stage variable, go to start of loop
      if(!isFirstStageIntFeas(primalSoln)) {
	if (0 == mype) cout << "Branching on first stage!\n";
	branchOnFirstStage(primalSoln);
	continue;
      }

      // If we get to this point, we know that the first stage decision variables
      // are integer feasible, but the LP solution is not integer feasible, so
      // one of the second stage scenarios must not be integer feasible, and
      // one of the variables in one of those scenarios should be branched on.
      if (0 == mype) cout << "Branching on second stage!\n";
      branchOnSecondStage(primalSoln);
      continue;
    }

    if (0 == mype) cout << "Objective function value = " << objUB << endl;
    if (0 == mype) cout << "Objective function LB = " << objLB << endl;

  }

  private:

  // Make default constructor impossible to call.
  BranchAndBoundTree();

} ;



int main(int argc, char **argv) {

        // Initialize MPI
	MPI_Init(&argc, &argv);

        // Get MPI process rank
	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

        // Help information if not enough arguments
	if (argc < 2) {
		if (0 == mype) printf("Usage: %s [SMPS root name]\n",argv[0]);
		return 1;
	}

	// Set PIPS logging level.
	PIPSLogging::init_logging(1);

        // Get SMPS file name and open SMPS file
	string smpsrootname(argv[1]);

	if (0 == mype) cout << "Reading SMPS input!\n";
	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

	//scoped_ptr<SMPSInput> s(new SMPSInput(datarootname,nscen));

        // Pass communicator to block angular data structures for data distribution
	BAContext ctx(MPI_COMM_WORLD);

	// Initialize branch-and-bound tree
	if (0 == mype) cout << "Initializing branch-and-bound tree!\n";
	BranchAndBoundTree bb(input);

	/*
	// Solve deterministic LP formulation via dual simplex
	PIPSSInterface solver(input, ctx, PIPSSInterface::useDual);

	if (argc == 5) {
		solver.loadStatus(argv[4]);
	}

	//solver.setDumpFrequency(5000,argv[3]);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();

	// Write solution (if given enough input arguments)
	if (argc >= 4 && argv[3][0] != '-') {
		if (0 == mype) printf("Writing solution\n");
		solver.writeStatus(argv[3]);
		if (0 == mype) printf("Finished writing solution\n");
	}
	*/

	if (bb.status != ProvenInfeasible) {
	  if (0 == mype) cout << "Calling branch-and-bound!\n";
	  bb.branchAndBound();
	}
	else { // bb.status == ProvenInfeasible
	  if (0 == mype) cout << "Problem is infeasible!\n";
	}


	// Clean up MPI data structures
	MPI_Finalize();

	return 0;
}

