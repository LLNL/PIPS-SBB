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
