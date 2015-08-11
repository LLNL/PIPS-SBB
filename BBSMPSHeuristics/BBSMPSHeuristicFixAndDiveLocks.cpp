#include "BBSMPSHeuristicFixAndDiveLocks.hpp"

using namespace std;


double objContribution(double valueToRound, double objCoefficient, int roundingDirection){
	
	double obj=fabs(objCoefficient);
	if (roundingDirection==0){//We are rounding down

		double differential= floor(valueToRound) - valueToRound;
		// cout<<"Obj contribution "<<valueToRound<< " O: "<<objCoefficient<<" RDir "<<roundingDirection<<" RESULT "<<objCoefficient*differential<<endl;
		
		return obj*differential;
	}
		
		double differential= ceil(valueToRound) - valueToRound;
	//	cout<<"Obj contribution "<<valueToRound<< " O: "<<objCoefficient<<" RDir "<<roundingDirection<<" RESULT "<<objCoefficient*differential<<endl;
		
		return obj*differential;
	


}

void generateLocks(denseBAVector & upLocks, denseBAVector &downLocks){

	
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	
	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
   
   	int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();
   	const BAFlagVector<constraintType> varTypes = rootSolver.getVariableTypes();
   	//cout<<"gl1"<<endl;
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int c = 0; c < input.nSecondStageCons(scen); c++)
			{

				const CoinShallowPackedVector row=rootSolver.retrieveTRow(c,scen);
				int nElems=row.getNumElements();
				const int*indices=row.getIndices();
				const double *elems=row.getElements(); 


		    	const CoinShallowPackedVector row2=rootSolver.retrieveWRow(c,scen);
		    	int nElems2=row2.getNumElements();
				const int* indices2=row2.getIndices();
				const double *elems2=row2.getElements(); 

				if (varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==LB){
					for (int el=0; el<nElems; el++){
			    		
			    		if (elems[el]<0)upLocks.getFirstStageVec()[indices[el]]++;
			    		else downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		if (elems2[el]<0)upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		else downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }
			    else if (varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==UB){
			    	for (int el=0; el<nElems; el++){
			    		
			    		if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
			    		else downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		if (elems2[el]>0)upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		else downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }
			    else if(varTypes.getSecondStageVec(scen)[input.nSecondStageVars(scen)+c]==Fixed){
			    	for (int el=0; el<nElems; el++){
			    		
			    		upLocks.getFirstStageVec()[indices[el]]++;
			    	    downLocks.getFirstStageVec()[indices[el]]++;

			    	}

			    	for (int el=0; el<nElems2; el++){
			    		upLocks.getSecondStageVec(scen)[indices2[el]]++;
			    		downLocks.getSecondStageVec(scen)[indices2[el]]++;

			    	}
			    }


		    }

		}
    }
    //cout<<"gl2"<<endl;

    //At this point each vector has its own count of first stage locks. Let's reduce
    double *upLock1stStagePtr = upLocks.getFirstStageVec().getPointer();
    double *downLock1stStagePtr = downLocks.getFirstStageVec().getPointer();
    MPI_Allreduce(MPI_IN_PLACE,upLock1stStagePtr,upLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,downLock1stStagePtr,downLocks.getFirstStageVec().length(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    //cout<<"gl3"<<endl;
    for (int c=0; c< firstStageRows ; c++){

    	const CoinShallowPackedVector row=rootSolver.retrieveARow(c);
    	int nElems=row.getNumElements();
		const int*indices=row.getIndices();
		const double *elems=row.getElements(); 

		if (varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==LB){
			for (int el=0; el<nElems; el++){
	    		if (elems[el]<0)upLocks.getFirstStageVec()[indices[el]]++;
	    		else downLocks.getFirstStageVec()[indices[el]]++;

	    	}
	    }
	    else if (varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==UB){
	    	for (int el=0; el<nElems; el++){
    			if (elems[el]>0)upLocks.getFirstStageVec()[indices[el]]++;
    			else downLocks.getFirstStageVec()[indices[el]]++;
    		}
	    }
	    else if(varTypes.getFirstStageVec()[input.nFirstStageVars()+c]==Fixed){
	    	for (int el=0; el<nElems; el++){
				upLocks.getFirstStageVec()[indices[el]]++;
	    	    downLocks.getFirstStageVec()[indices[el]]++;

	    	}
    	}

    }
    //cout<<"gl4"<<endl;
}



bool BBSMPSHeuristicFixAndDiveLocks::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB){
	
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Fix and Dive heuristic.";
	timesCalled++;


	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	
	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const denseBAVector varObjectives = rootSolver.getVarObjective();
    

    const BAFlagVector<constraintType> varTypes = rootSolver.getVariableTypes();

	
    int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();
   
   	int MAX_ITERS=firstStageVars;
   	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		MAX_ITERS+=input.nSecondStageVars(scen);
	}


	denseBAVector upLocks;
    denseBAVector downLocks;
    upLocks.allocate(dimsSlacks, ctx, PrimalVector);
	downLocks.allocate(dimsSlacks, ctx, PrimalVector);
   	upLocks.clear();
   	downLocks.clear();
    
    generateLocks(upLocks,downLocks);

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
	double objCont=COIN_DBL_MAX;//We are trying to minimize
			
	//Order by most to least fractional
    for (int i=0; i< input.nFirstStageVars(); i++){
    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=upLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];
    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0);
    		
    	}
    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
    		bestLockIndex=i;
    		bestLock=downLocks.getFirstStageVec()[i];
    		bestVarObj=variableObjectives.getFirstStageVec()[i];
    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1);
    	
    	}
    }
    bool allInteger=(bestLockIndex==-1);
	int iter=0;
	while(!allInteger && iter<MAX_ITERS){
		
		if (bestLock==upLocks.getFirstStageVec()[bestLockIndex] && bestLock==downLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=roundToNearestInteger(auxSolution.getFirstStageVec()[bestLockIndex]);
		//	cout<<"Rounding Nearest -1 "<<bestLockIndex<<" to "<<lb.getFirstStageVec()[bestLockIndex]<<" Score: "<<downLocks.getFirstStageVec()[bestLockIndex]<<" "<<upLocks.getFirstStageVec()[bestLockIndex]<<" obj "<<variableObjectives.getFirstStageVec()[bestLockIndex]<<endl;
	
		}
		else if (bestLock==upLocks.getFirstStageVec()[bestLockIndex]){

			lb.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=floor(auxSolution.getFirstStageVec()[bestLockIndex]);	
		//	cout<<"Rounding down -1 "<<bestLockIndex<<" to "<<lb.getFirstStageVec()[bestLockIndex]<<" Score: "<<downLocks.getFirstStageVec()[bestLockIndex]<<" "<<upLocks.getFirstStageVec()[bestLockIndex]<<" obj "<<variableObjectives.getFirstStageVec()[bestLockIndex]<<endl;
		}
		else{

			lb.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=ceil(auxSolution.getFirstStageVec()[bestLockIndex]);
		//	cout<<"Rounding up -1 "<<bestLockIndex<<" to "<<lb.getFirstStageVec()[bestLockIndex]<<" Score: "<<downLocks.getFirstStageVec()[bestLockIndex]<<" "<<upLocks.getFirstStageVec()[bestLockIndex]<<" obj "<<variableObjectives.getFirstStageVec()[bestLockIndex]<<endl;
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
	    objCont=COIN_DBL_MAX;//We are trying to minimize
			
		//Order by most to least fractional
	    for (int i=0; i< input.nFirstStageVars(); i++){
	    	if (((upLocks.getFirstStageVec()[i]> bestLock)|| (upLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestLockIndex=i;
	    		bestLock=upLocks.getFirstStageVec()[i];
	    		bestVarObj=variableObjectives.getFirstStageVec()[i];
	    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],0);
	    		
	    	}
	    	if (((downLocks.getFirstStageVec()[i]> bestLock)|| (downLocks.getFirstStageVec()[i]== bestLock && objCont> objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1))) && input.isFirstStageColInteger(i) && !isIntFeas(auxSolution.getFirstStageVec()[i], intTol)){
	    		bestLockIndex=i;
	    		bestLock=downLocks.getFirstStageVec()[i];
	    		bestVarObj=variableObjectives.getFirstStageVec()[i];
	    		objCont= objContribution(auxSolution.getFirstStageVec()[i],varObjectives.getFirstStageVec()[i],1);
	    	
	    	}
	    }
	    iter++;
	    allInteger=(bestLockIndex==-1);

	}
 

	int bestLockScen=-1;
	bestLockIndex=-1;
	bestLock=-1;
	bestVarObj=COIN_DBL_MAX;
	 objCont=COIN_DBL_MAX;//We are trying to minimize
	
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int i = 0; i < input.nSecondStageVars(scen); i++)
			{
				if (((upLocks.getSecondStageVec(scen)[i] >bestLock) || (upLocks.getSecondStageVec(scen)[i] ==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=upLocks.getSecondStageVec(scen)[i];
		    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0);
		    	}
		    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
		    		bestLockIndex=i;
		    		bestLockScen=scen;
		    		bestLock=downLocks.getSecondStageVec(scen)[i];
		    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1);
		    	}
		    }
		}
	}
	
	
	int maxCont;
	int errorFlag = MPI_Allreduce(&bestLockIndex, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());
	int iteration=0;
	iter =0;
	while (maxCont>-1 && iteration<MAX_ITERS){
		iteration++;
		if (bestLockIndex>-1){
			// cout<<input.nSecondStageVars(bestFracScen)<<" we are about to update!!"<<mype<<" has chosen "<<bestFracIndex<<" frac part "<<bestFracPart<<" frac scen "<<bestFracScen<<" "<<maxCont<<endl;
	
		if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] && bestLock==downLocks.getSecondStageVec(bestLockScen)[bestLockIndex]){
			lb.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			ub.getSecondStageVec(bestLockScen)[bestLockIndex]=roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
		//	 cout<<" we are about to update!!"<<mype<<" has chosen "<<bestLockIndex<<" frac  "<<roundToNearestInteger(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex])<<"best lock "<<bestLock<<endl;
	
		}

	
			if (bestLock==upLocks.getSecondStageVec(bestLockScen)[bestLockIndex] ){
		//		 cout<<" we are about to update!!"<<mype<<" has chosen "<<bestLockIndex<<" frac  "<<floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex])<<"best lock "<<bestLock<<endl;
	
				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			
			}
			else{
		//		 cout<<" we are about to update!!"<<mype<<" has chosen "<<bestLockIndex<<" frac  "<<ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex])<<"best lock "<<bestLock<<endl;
	
				lb.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
				ub.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(auxSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			}
		
			
		}
		rootSolver.setLB(lb);
			rootSolver.setUB(ub);
		//cout<<mype<<" got here "<<endl;
		rootSolver.commitStates();	
//cout<<mype<<" got here 2"<<endl;
		rootSolver.go();
//cout<<mype<<" got here3 "<<endl;
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
		objCont=COIN_DBL_MAX;//We are trying to minimize
		
		for (int scen = 0; scen < input.nScenarios(); scen++)
		{
			if(ctx.assignedScenario(scen)) {
				for (int i = 0; i < input.nSecondStageVars(scen); i++)
				{
					if (((upLocks.getSecondStageVec(scen)[i] >bestLock) || (upLocks.getSecondStageVec(scen)[i] ==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		bestLockIndex=i;
			    		bestLockScen=scen;
			    		bestLock=upLocks.getSecondStageVec(scen)[i];
			    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],0);
			    	}
			    	if (((downLocks.getSecondStageVec(scen)[i] >bestLock)||(downLocks.getSecondStageVec(scen)[i]==bestLock && objCont> objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1))) && input.isSecondStageColInteger(scen,i) && !isIntFeas(auxSolution.getSecondStageVec(scen)[i], intTol)){
			    		bestLockIndex=i;
			    		bestLockScen=scen;
			    		bestLock=downLocks.getSecondStageVec(scen)[i];
			    		objCont= objContribution(auxSolution.getSecondStageVec(scen)[i],varObjectives.getSecondStageVec(scen)[i],1);
			    	}
			    }
			}
		}
		
		maxCont=-2;
	    errorFlag = MPI_Allreduce(&bestLock, &maxCont, 1, MPI_INT,  MPI_MAX, ctx.comm());
	    iter++;
	   // cout<<iteration<<" "<<mype<<" has chosen "<<bestLockIndex<<" frac part "<<bestLock<<" frac scen "<<bestLockScen<<" "<<maxCont<<endl;

	}

   	
	solverState lpStatus = rootSolver.getStatus();
	bool otherThanOptimal = (Optimal != lpStatus); 
	
	if(!otherThanOptimal){
		denseBAVector solVector=rootSolver.getPrimalSolution();
		solution=BBSMPSSolution(solVector,rootSolver.getObjective());

		cout<<"DID WE FIND A LOCK ROUNDING SOLUTION "<<isLPIntFeas(solVector)<<" of qual "<<rootSolver.getObjective()<<endl;

		 
	}
	//return if success
	bool success= (!otherThanOptimal && rootSolver.getObjective()<objUB);
	timesSuccessful+=(success);
	
	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicFixAndDiveLocks::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	if (node->getNodeDepth()<10)return true;
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	
	int numberOfFractionalVariables=0;
	int nIntVars=0;
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			if (!isIntFeas(LPRelaxationSolution.getFirstStageVec()[col],intTol)) numberOfFractionalVariables++;
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
					if (!isIntFeas(LPRelaxationSolution.getSecondStageVec(scen)[col],intTol)) numberOfFractionalVariables2++;
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
	cout<<"Total number of frac vars "<<numberOfFractionalVariables<<" n int vars "<<nIntVars<<	" ratio "<<(numberOfFractionalVariables*100/nIntVars)<<endl;
	return ((numberOfFractionalVariables*100/nIntVars)<40 );

}