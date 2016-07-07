#include "BBSMPSHeuristicFixAndDive.hpp"

using namespace std;


bool variableVectorSort(pair <int, double> i,pair <int, double> j) { return (i.second<j.second);}
bool secondStageVariableVectorSort(pair <pair<int, int> , double> i,pair <pair<int, int> , double> j) { return (i.second<j.second);}



bool BBSMPSHeuristicFixAndDive::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();

	timesCalled++;
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Fix and Dive heuristic.";


	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());

	node->getAllBranchingInformation(lb,ub);

	denseBAVector auxSolution(LPRelaxationSolution);
	//Get all variables
    double bestFracPart=2;
    int bestFracIndex=-1;


	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	node->reconstructWarmStartState(ps);
		rootSolver.setStates(ps);
	//Order by most to least fractional
    for (int i=0; i< input.nFirstStageVars(); i++){
    	if (fracPart(auxSolution.getFirstStageVec()[i])<bestFracPart && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestFracIndex=i;
    		bestFracPart=fracPart(auxSolution.getFirstStageVec()[i]);
    	}
    }
    bool allInteger=(bestFracIndex==-1);

	while(!allInteger){
		lb.getFirstStageVec()[bestFracIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestFracIndex]);
		ub.getFirstStageVec()[bestFracIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestFracIndex]);

		rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		rootSolver.commitStates();

		rootSolver.go();

		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus);
		if (otherThanOptimal) {
			cumulativeTime+=(MPI_Wtime()-startTimeStamp);
			return false;
		}
		auxSolution=rootSolver.getPrimalSolution();


		bestFracPart=2;
     	bestFracIndex=-1;

    	//Order by most to least fractional
	    for (int i=0; i< input.nFirstStageVars(); i++){
	    	if (fracPart(auxSolution.getFirstStageVec()[i])<bestFracPart && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestFracIndex=i;
	    		bestFracPart=fracPart(auxSolution.getFirstStageVec()[i]);
	    	}
	    }
	    allInteger=(bestFracIndex==-1);

	}

	for (int i=0; i< input.nFirstStageVars(); i++){
 		lb.getFirstStageVec()[i]=auxSolution.getFirstStageVec()[i];
		ub.getFirstStageVec()[i]=auxSolution.getFirstStageVec()[i];
	}

	int bestFracScen=-1;
	bestFracIndex=-1;
	bestFracPart=2;
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int i = 0; i < input.nSecondStageVars(scen); i++)
			{
				if (fracPart(auxSolution.getSecondStageVec(scen)[i])<bestFracPart && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestFracIndex=i;
		    		bestFracScen=scen;
		    		bestFracPart=fracPart(auxSolution.getSecondStageVec(scen)[i]);
		    	}
		    }
		}
	}


	int maxCont;
	int errorFlag = MPI_Allreduce(&bestFracIndex, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());
	int iteration=0;
	while (maxCont>-1){
		iteration++;
		if (bestFracIndex>-1){

			lb.getSecondStageVec(bestFracScen)[bestFracIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestFracScen)[bestFracIndex]);
			ub.getSecondStageVec(bestFracScen)[bestFracIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestFracScen)[bestFracIndex]);


		}
		rootSolver.setLB(lb);
		rootSolver.setUB(ub);
		rootSolver.commitStates();
		rootSolver.go();
		solverState lpStatus = rootSolver.getStatus();
		int otherThanOptimal = (Optimal != lpStatus);
		int anyOtherThanOptimal=0;
	     MPI_Allreduce(&otherThanOptimal, &anyOtherThanOptimal, 1, MPI_INT,  MPI_MAX, ctx.comm());

		if (anyOtherThanOptimal!=0){
			cumulativeTime+=(MPI_Wtime()-startTimeStamp);
			return false;
		}
		auxSolution=rootSolver.getPrimalSolution();


		bestFracPart=2;
     	bestFracIndex=-1;

    	//Order by most to least fractional
	    bestFracScen=-1;
		bestFracIndex=-1;
		bestFracPart=2;

		int fracVals=0;

		for (int scen = 0; scen < input.nScenarios(); scen++)
		{
			if(ctx.assignedScenario(scen)) {
				for (int i = 0; i < input.nSecondStageVars(scen); i++)
				{
					if (fracPart(auxSolution.getSecondStageVec(scen)[i])<bestFracPart && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		bestFracIndex=i;
			    		bestFracScen=scen;
			    		bestFracPart=fracPart(auxSolution.getSecondStageVec(scen)[i]);
			    	}
			    	if (input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		fracVals++;
			    	}
			    }
			}
		}


		maxCont=0;
	    errorFlag = MPI_Allreduce(&bestFracIndex, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());

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

bool BBSMPSHeuristicFixAndDive::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
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
	return ((numberOfFractionalVariables*100/nIntVars)<25 );

}