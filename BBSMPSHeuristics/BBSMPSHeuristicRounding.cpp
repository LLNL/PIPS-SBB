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
#include "BBSMPSHeuristicRounding.hpp"

using namespace std;


bool BBSMPSHeuristicRounding::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	timesCalled++;
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the simple rounding heuristic.";


	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	node->getAllBranchingInformation(lb,ub);

	//Apply Simple rounding to 1st stage vars
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{
		if(input.isFirstStageColInteger(col)){
			lb.getFirstStageVec()[col]=roundToNearestInteger(LPRelaxationSolution.getFirstStageVec()[col]);
			ub.getFirstStageVec()[col]=roundToNearestInteger(LPRelaxationSolution.getFirstStageVec()[col]);
		}
	}

	rootSolver.setLB(lb);
	rootSolver.setUB(ub);

	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	node->reconstructWarmStartState(ps);
	rootSolver.setStates(ps);
	rootSolver.commitStates();


	//Check if feasible
	rootSolver.go();

	/* Check solver status for infeasibility/optimality */
	solverState lpStatus = rootSolver.getStatus();
	bool otherThanOptimal = (Optimal != lpStatus);
	if (otherThanOptimal) return false;

	denseBAVector primalSoln(rootSolver.getPrimalSolution());
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	//Apply simple rounding to 2nd stage vars
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				if(input.isSecondStageColInteger(scen,col)){
					lb.getSecondStageVec(scen)[col]=roundToNearestInteger(primalSoln.getSecondStageVec(scen)[col]);
					ub.getSecondStageVec(scen)[col]=roundToNearestInteger(primalSoln.getSecondStageVec(scen)[col]);
				}
			}
		}
	}

	//Check if feasible
	rootSolver.setLB(lb);
	rootSolver.setUB(ub);
	//Check if feasible
	rootSolver.setStates(ps);
	rootSolver.commitStates();
	rootSolver.go();
	lpStatus = rootSolver.getStatus();
	otherThanOptimal = (Optimal != lpStatus);

	double objUB=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();


	if(!otherThanOptimal && rootSolver.getObjective()<objUB){
		denseBAVector solVector=rootSolver.getPrimalSolution();
		if (isLPIntFeas(solVector)){
			BBSMPSSolution sol(solVector,rootSolver.getObjective());
			sol.setTimeOfDiscovery(BBSMPSSolver::instance()->getWallTime());
			BBSMPSSolver::instance()->addSolutionToPool(sol);
		}

	}
	//return if success
	bool success= (!otherThanOptimal && rootSolver.getObjective()<objUB);
	timesSuccessful+=success;

	if (0 == mype && success) BBSMPS_ALG_LOG_SEV(info) << "The simple rounding heuristic was successful.";

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicRounding::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();

	int numberOfFractionalVariables=0;
	int nIntVars=0;
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{
		if(input.isFirstStageColInteger(col)){
			numberOfFractionalVariables+=(!isIntFeas(LPRelaxationSolution.getFirstStageVec()[col],intTol));
			nIntVars++;
		}

	}

	int numberOfFractionalVariables2=0;
	int nIntVars2=0;


	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{

				if(input.isSecondStageColInteger(scen,col)){
					numberOfFractionalVariables2+=(!isIntFeas(LPRelaxationSolution.getSecondStageVec(scen)[col],intTol));
					nIntVars2++;

				}
			}
		}
	}

	int totalCount2;
	int errorFlag = MPI_Allreduce(&numberOfFractionalVariables2,
		&totalCount2,
		1,
		MPI_INT,
		MPI_SUM,
		ctx.comm());

	int totalIntVars2;
	errorFlag = MPI_Allreduce(&nIntVars2,
		&totalIntVars2,
		1,
		MPI_INT,
		MPI_SUM,
		ctx.comm());


	nIntVars+=totalIntVars2;
	numberOfFractionalVariables=+totalCount2;
	if (nIntVars==0)return false;
	return ((numberOfFractionalVariables*100/nIntVars)<10 );

}
