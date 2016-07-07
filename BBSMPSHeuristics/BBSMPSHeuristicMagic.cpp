#include "BBSMPSHeuristicMagic.hpp"

using namespace std;




bool BBSMPSHeuristicMagic::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	int mype=BBSMPSSolver::instance()->getMype();

	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Magic heuristic.";

	double startTimeStamp = MPI_Wtime();

	timesCalled++;

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();



    int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();


	denseBAVector upLocks;
    denseBAVector downLocks;
    upLocks.allocate(dimsSlacks, ctx, PrimalVector);
	downLocks.allocate(dimsSlacks, ctx, PrimalVector);
   	upLocks.clear();
   	downLocks.clear();

    for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int c = 0; c < input.nSecondStageCons(scen); c++)
			{
				const CoinShallowPackedVector row=rootSolver.retrieveTRow(c,scen);
				int nElems=row.getNumElements();
				const int*indices=row.getIndices();
				const double *elems=row.getElements();
		    	for (int el=0; el<nElems; el++){
		    		if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
		    		else downLocks.getFirstStageVec()[indices[el]]++;

		    	}

		    	const CoinShallowPackedVector row2=rootSolver.retrieveWRow(c,scen);
		    	int nElems2=row2.getNumElements();
				const int* indices2=row2.getIndices();
				const double *elems2=row2.getElements();

		    	for (int el=0; el<nElems2; el++){
		    		if (elems2[el]>0)upLocks.getSecondStageVec(scen)[indices2[el]]++;
		    		else downLocks.getSecondStageVec(scen)[indices2[el]]++;

		    	}


		    }

		}
    }


    //At this point each vector has its own count of first stage locks. Let's reduce
    double *upLock1stStagePtr = upLocks.getFirstStageVec().getPointer();
    double *downLock1stStagePtr = downLocks.getFirstStageVec().getPointer();
    MPI_Allreduce(MPI_IN_PLACE,upLock1stStagePtr,upLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,downLock1stStagePtr,downLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


    for (int c=0; c< firstStageRows ; c++){

    	const CoinShallowPackedVector row=rootSolver.retrieveARow(c);
    	int nElems=row.getNumElements();
		const int*indices=row.getIndices();
		const double *elems=row.getElements();

    	for (int el=0; el<nElems; el++){
    		if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
    		else downLocks.getFirstStageVec()[indices[el]]++;

    	}


    }

    denseBAVector variableObjectives = rootSolver.getVarObjective();
   	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());

	node->getAllBranchingInformation(lb,ub);

	denseBAVector auxSolution(LPRelaxationSolution);

	//Get all variables
	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	node->reconstructWarmStartState(ps);
	rootSolver.setStates(ps);

	double bestLockIndex=-1;
    int bestLock=-1;
    double bestVarObj=COIN_DBL_MAX;
	//Order by most to least fractional
    for (int i=0; i< input.nFirstStageVars(); i++){
    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && bestVarObj>variableObjectives.getFirstStageVec()[i])) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=upLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];

    	}
    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && bestVarObj>variableObjectives.getFirstStageVec()[i])) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=downLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];

    	}
    }
    bool allInteger=(bestLockIndex==-1);

	while(!allInteger){

		if (bestLock==upLocks.getFirstStageVec()[bestLockIndex] && bestLock==downLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);

		}
		else if (bestLock==upLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);
		}
		else{

			lb.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
		}

		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		rootSolver.commitStates();

		rootSolver.go();

		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus);
		if (otherThanOptimal) return false;
		auxSolution=rootSolver.getPrimalSolution();

		bestLockIndex=-1;
  		bestLock=-1;
		bestVarObj=COIN_DBL_MAX;
	    for (int i=0; i< input.nFirstStageVars(); i++){
	    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && bestVarObj>variableObjectives.getFirstStageVec()[i])) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=upLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];

	    	}
	    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && bestVarObj>variableObjectives.getFirstStageVec()[i])) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestLockIndex=i;
	    		bestLock=downLocks.getFirstStageVec()[i];
	    		bestVarObj=variableObjectives.getFirstStageVec()[i];

	    	}
	    }
	    allInteger=(bestLockIndex==-1);

	}


	int bestLockScen=-1;
	bestLockIndex=-1;
	bestLock=-1;
	bestVarObj=COIN_DBL_MAX;
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int i = 0; i < input.nSecondStageVars(scen); i++)
			{
				if (((upLocks.getSecondStageVec(scen)[i] >bestLock) || (upLocks.getSecondStageVec(scen)[i] ==bestLock && bestVarObj>variableObjectives.getSecondStageVec(scen)[i])) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=upLocks.getSecondStageVec(scen)[i];
		    	}
		    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && bestVarObj>variableObjectives.getSecondStageVec(scen)[i])) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=downLocks.getSecondStageVec(scen)[i];
		    	}
		    }
		}
	}


	int maxCont;
	int errorFlag = MPI_Allreduce(&bestLockIndex, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());
	int iteration=0;
	while (maxCont>-1){
		iteration++;
		if (bestLockIndex>-1){

		if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] && bestLock==downLocks.getSecondStageVec(bestLockScen)[bestLockIndex]){
			lb.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			ub.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);

		}


			if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] ){

				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);

			}
			else{

				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			}


		}
		rootSolver.setLB(lb);
			rootSolver.setUB(ub);
		rootSolver.commitStates();
		rootSolver.go();
		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus);
		int anyOtherThanOptimal=0;
	     MPI_Allreduce(&anyOtherThanOptimal, &otherThanOptimal, 1, MPI_INT,  MPI_MAX, ctx.comm());

		if (anyOtherThanOptimal>0) return false;
		auxSolution=rootSolver.getPrimalSolution();


		bestLockScen=-1;
		bestLockIndex=-1;
		bestLock=-1;

		bestVarObj=COIN_DBL_MAX;
		for (int scen = 0; scen < input.nScenarios(); scen++)
		{
			if(ctx.assignedScenario(scen)) {
				for (int i = 0; i < input.nSecondStageVars(scen); i++)
				{
					if (((upLocks.getSecondStageVec(scen)[i] >bestLock)||(upLocks.getSecondStageVec(scen)[i]==bestLock && bestVarObj>variableObjectives.getSecondStageVec(scen)[i])) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		bestLockIndex=i;
			    		bestLockScen=scen;
			    		bestLock=upLocks.getSecondStageVec(scen)[i];
			    	}
			    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && bestVarObj>variableObjectives.getSecondStageVec(scen)[i])) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		bestLockIndex=i;
			    		bestLockScen=scen;
			    		bestLock=downLocks.getSecondStageVec(scen)[i];
			    	}
			    }
			}
		}

		maxCont=-2;
	    errorFlag = MPI_Allreduce(&bestLock, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());


	}


	solverState lpStatus = rootSolver.getStatus();
	bool otherThanOptimal = (Optimal != lpStatus);

	double objUB=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	if(!otherThanOptimal){
		denseBAVector solVector=rootSolver.getPrimalSolution();
		if (isLPIntFeas(solVector)){
					BBSMPSSolution sol(solVector,rootSolver.getObjective());
					sol.setTimeOfDiscovery(BBSMPSSolver::instance()->getWallTime());
					BBSMPSSolver::instance()->addSolutionToPool(sol);
		}


	}
	//return if success
	bool success= (!otherThanOptimal && rootSolver.getObjective()<objUB);
	timesSuccessful+=(success);


	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicMagic::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	/*SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
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
	return ((numberOfFractionalVariables*100/nIntVars)<25 );
	*/
	if (node->getNodeDepth()<10)return true;
	int nodeDepth=node->getNodeDepth();
	assert(nodeDepth>=0);
	return (nodeDepth%15==0);
}