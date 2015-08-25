// ----------------------------------------------------------------------------
/**
   File: BBSMPSSolverState.hpp

   Description: The solver state class hides the complexity behind managing all the
   				possible solver states.

   Limitations: All logic regarding the state of the solver should be allocated only
   				in this class.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSSOLVERSTATE_H
#define BBSMPSSOLVERSTATE_H
#include "PIPSSInterface.hpp"


class BBSMPSSolverState {

public:
	BBSMPSSolverState(solverState s);
	
	~BBSMPSSolverState();

	void setStatusToPrimalFeasible() ;

	void setStatusToOptimal() ;

	void setStatusToProvenInfeasible() ;

	void setStatusToStopped() ;

	bool isPrimalFeasible();

	bool isLoadedFromFile();

private:
	BBSMPSSolverState();

	solverState status;

};

#endif