#include "BBSMPSSolverState.hpp"

BBSMPSSolverState::BBSMPSSolverState(solverState s){
	status=s;
}

BBSMPSSolverState::~BBSMPSSolverState() {}
void BBSMPSSolverState::setStatusToPrimalFeasible() {
	bool isInReachableState = (LoadedFromFile == status);
		    //      || (status == Bounded);
	if (isInReachableState) {
		status = PrimalFeasible;
	}
}

void BBSMPSSolverState::setStatusToOptimal() {
	bool isInReachableState = (LoadedFromFile == status) || (PrimalFeasible == status); // || (Bounded == status);
	if (isInReachableState) {
		status = Optimal;
	}
}

void BBSMPSSolverState::setStatusToProvenInfeasible() {
	bool isInReachableState = (LoadedFromFile == status) ||
	(ProvenInfeasible == status);
	if (isInReachableState) {
		status = ProvenInfeasible;
	}
}

void BBSMPSSolverState::setStatusToStopped() {
	bool isInReachableState = (LoadedFromFile == status);
	if (isInReachableState) {
		status = Stopped;
	}
}


bool BBSMPSSolverState::isPrimalFeasible(){
	return (status==PrimalFeasible);
}

bool BBSMPSSolverState::isLoadedFromFile(){
	return (status==LoadedFromFile);
}
