#include "BBSMPSHeuristicLockRounding.hpp"

using namespace std;


double objContribution2(double valueToRound, double objCoefficient, int roundingDirection){
	
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

void generateLocks2(denseBAVector & upLocks, denseBAVector &downLocks){

	
	
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
int findAndFixFirstStageConstraint(denseBAVector &roundedSolution,denseBAVector & upLocks, denseBAVector &downLocks,denseBAVector &lb,denseBAVector &ub){

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	int mype=BBSMPSSolver::instance()->getMype();
	
	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
   	const denseBAVector varObjectives = rootSolver.getVarObjective();
    
    int firstStageVars=input.nFirstStageVars();
    int firstStageRows=input.nFirstStageCons();

    //cout<<"ff1"<<endl;
	//First check first stage
	bool allRowsFeasible=true;
	

	for (int c=0; c< firstStageRows &&allRowsFeasible; c++){
		int brokenDirection=rootSolver.isRowFeasible(c, -1, roundedSolution);
		//cout<<"Trying to fix constraint "<<c<<" which is broken on the side "<<brokenDirection<<endl;
    		
		if (brokenDirection!=1){
			allRowsFeasible=false;
			const CoinShallowPackedVector row=rootSolver.retrieveARow(c);
	    	int nElems=row.getNumElements();
			const int*indices=row.getIndices();
			const double *elems=row.getElements(); 
			double bestLock=COIN_DBL_MAX;
			int varIndexOfBestLock=-1;
			double objCont=COIN_DBL_MAX;//We are trying to minimize
			int directionToRound;
			for (int el=0; el<nElems; el++){
	    		if (input.isFirstStageColInteger(indices[el]) && !isIntFeas(roundedSolution.getFirstStageVec()[indices[el]],intTol)){ //We have a possible rounding candidate
	    			if ((elems[el]>0 && brokenDirection==-1) || (elems[el]<0 && brokenDirection==-2) ){ ///If we have the constraint broken in the ub and the cofficient is positive, OR we have the constraint broken in the lb and the coefficient is negative we must round down
	    				if (bestLock>downLocks.getFirstStageVec()[indices[el]] || (bestLock==downLocks.getFirstStageVec()[indices[el]] && objCont> objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],0))){
	    					bestLock=downLocks.getFirstStageVec()[indices[el]];
	    					varIndexOfBestLock=indices[el];
	    					directionToRound=0;
	    					objCont=objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],0);
	    				}
	    			}
	    			else if ((elems[el]<0 && brokenDirection==-1) || (elems[el]>0 && brokenDirection==-2)){ ///If we have the constraint broken in the ub and the cofficient is negative, OR we have the constraint broken in the ub and the coefficient is positive we must round down
	    				if (bestLock>upLocks.getFirstStageVec()[indices[el]] || (bestLock==upLocks.getFirstStageVec()[indices[el]] && objCont> objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],1))){
	    					bestLock=upLocks.getFirstStageVec()[indices[el]];
	    					varIndexOfBestLock=indices[el];
	    					directionToRound=1;
	    					objCont=objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],1);
	    				}
	    			}
	    		}
	    		

	    	}
	    	if (varIndexOfBestLock != -1 && directionToRound==0){//Round down variable
	    		//cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock])<<endl;
	    		roundedSolution.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		lb.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		
	    		ub.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		
	    		return 2;
	    	}
	    	else if (varIndexOfBestLock != -1 && directionToRound==1){//Round up variable
	    		//cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock])<<endl;
	    		roundedSolution.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		lb.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		ub.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
	    		return 2;
	    	}

		}

	}

/*
	//Al first stage rows seem fine. Now we should look at the T rows of each scenario
	double bestLock=COIN_DBL_MAX;
	int varIndexOfBestLock=-1;
	int scenOfBestLock=-1;
	double objCont=COIN_DBL_MAX;//We are trying to minimize
	int directionToRound;

	for (int scen = 0; scen < input.nScenarios() && allRowsFeasible; scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int c = 0; c < input.nSecondStageCons(scen) && allRowsFeasible; c++)
			{
				int brokenDirection=rootSolver.isRowFeasible(c, scen, roundedSolution);
				////cout<<"Trying to fix constraint "<<c<<" s"<<scen<<" which is broken on the side "<<brokenDirection<<endl;
		    		
				if (brokenDirection!=1){
					allRowsFeasible=false;
					
					const CoinShallowPackedVector row=rootSolver.retrieveTRow(c,scen);					
			
			    	int nElems=row.getNumElements();
					const int*indices=row.getIndices();
					const double *elems=row.getElements(); 
				
					for (int el=0; el<nElems; el++){
			    		if (input.isFirstStageColInteger(indices[el]) && !isIntFeas(roundedSolution.getFirstStageVec()[indices[el]],intTol)){ //We have a possible rounding candidate
			    			if ((elems[el]>0 && brokenDirection==-1) || (elems[el]<0 && brokenDirection==-2) ){ ///If we have the constraint broken in the ub and the cofficient is positive, OR we have the constraint broken in the lb and the coefficient is negative we must round down
			    				if (bestLock>downLocks.getFirstStageVec()[indices[el]] || (bestLock==downLocks.getFirstStageVec()[indices[el]] && objCont> objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],0))){
			    					bestLock=downLocks.getFirstStageVec()[indices[el]];
			    					varIndexOfBestLock=indices[el];
			    					directionToRound=0;
			    					objCont=objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],0);
			    				}
			    			}
			    			else if ((elems[el]<0 && brokenDirection==-1) || (elems[el]>0 && brokenDirection==-2)){ ///If we have the constraint broken in the ub and the cofficient is negative, OR we have the constraint broken in the ub and the coefficient is positive we must round down
			    				if (bestLock>upLocks.getFirstStageVec()[indices[el]] || (bestLock==upLocks.getFirstStageVec()[indices[el]] && objCont> objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],1))){
			    					bestLock=upLocks.getFirstStageVec()[indices[el]];
			    					varIndexOfBestLock=indices[el];
			    					directionToRound=1;
			    					objCont=objContribution2(roundedSolution.getFirstStageVec()[indices[el]],varObjectives.getFirstStageVec()[indices[el]],1);
			    				}
			    			}
			    		}
			    		

			    	}
			    }
			}
		}
	}
	int my[2];
	int best[2];
	my[0]=bestLock;
	my[1]=mype;
	
	MPI_Allreduce(&my,&best,1,MPI_2INT,MPI_MAXLOC,ctx.comm());

	if (best[1]==mype && best[0]!=COIN_DBL_MAX){
		if (varIndexOfBestLock != -1 && directionToRound==0){//Round down variable
			////cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock])<<endl;
			roundedSolution.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			lb.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			
			ub.getFirstStageVec()[varIndexOfBestLock]=floor(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			
			return 2;
		}
		else if (varIndexOfBestLock != -1 && directionToRound==1){//Round up variable
			////cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock])<<endl;
			roundedSolution.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			lb.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			ub.getFirstStageVec()[varIndexOfBestLock]=ceil(roundedSolution.getFirstStageVec()[varIndexOfBestLock]);
			return 2;
		}

	}
	else if (best[0]!=COIN_DBL_MAX) return 2;
	
	
*/

	return allRowsFeasible;
}

int findFreshFirstStageVar(denseBAVector &roundedSolution,denseBAVector & upLocks, denseBAVector &downLocks,denseBAVector &lb,denseBAVector &ub){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	int bestLockIndex=-1;
	int bestLockScen=-1;
	double maxLock=COIN_DBL_MIN;
	double objCont=COIN_DBL_MAX;//We are trying to minimize
		PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
   		
	const denseBAVector varObjectives = rootSolver.getVarObjective();
    //cout<<"FF1"<<endl;
	int roundDirection;
	for (int v=0; v<input.nFirstStageVars(); v++){
		if(input.isFirstStageColInteger(v) && !isIntFeas(roundedSolution.getFirstStageVec()[v],intTol)){
			//We have a possible candidate
			if(maxLock<upLocks.getFirstStageVec()[v] || (maxLock==upLocks.getFirstStageVec()[v] && objCont> objContribution2(roundedSolution.getFirstStageVec()[v],varObjectives.getFirstStageVec()[v],0) )){
				maxLock=upLocks.getFirstStageVec()[v];
				bestLockIndex=v;
				roundDirection=0;
				objCont=objContribution2(roundedSolution.getFirstStageVec()[v],varObjectives.getFirstStageVec()[v],0);
	    				
			}
			if(maxLock<downLocks.getFirstStageVec()[v]|| (maxLock==downLocks.getFirstStageVec()[v] && objCont> objContribution2(roundedSolution.getFirstStageVec()[v],varObjectives.getFirstStageVec()[v],1) )){
				maxLock=downLocks.getFirstStageVec()[v];
				bestLockIndex=v;
				roundDirection=1;
				objCont=objContribution2(roundedSolution.getFirstStageVec()[v],varObjectives.getFirstStageVec()[v],1);
	    		
			}
		}
	}

	if (bestLockIndex!=-1){
		if (roundDirection==1) {
			//cout<<"rounded var "<<bestLockIndex<<" "<<maxLock<<" to be "<<ceil(roundedSolution.getFirstStageVec()[bestLockIndex])<<endl;	
			roundedSolution.getFirstStageVec()[bestLockIndex]=ceil(roundedSolution.getFirstStageVec()[bestLockIndex]);
			lb.getFirstStageVec()[bestLockIndex]=ceil(roundedSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=ceil(roundedSolution.getFirstStageVec()[bestLockIndex]);
			return 1;
		}
		else{
			//cout<<"rounded var "<<bestLockIndex<<" "<<maxLock<<" to be "<<floor(roundedSolution.getFirstStageVec()[bestLockIndex])<<endl;
			roundedSolution.getFirstStageVec()[bestLockIndex]=floor(roundedSolution.getFirstStageVec()[bestLockIndex]);
			lb.getFirstStageVec()[bestLockIndex]=floor(roundedSolution.getFirstStageVec()[bestLockIndex]);
			ub.getFirstStageVec()[bestLockIndex]=floor(roundedSolution.getFirstStageVec()[bestLockIndex]);
			return 1;
		} 
	}
	return 0;
}
int findAndFixSecondStageConstraint(denseBAVector &roundedSolution,denseBAVector & upLocks, denseBAVector &downLocks,denseBAVector &lb,denseBAVector &ub){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	
	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const denseBAVector varObjectives = rootSolver.getVarObjective();
   
	//First check first stage
	bool allRowsFeasible=true;
	//cout<<"ff2"<<endl;

	for (int scen = 0; scen < input.nScenarios() && allRowsFeasible; scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int c = 0; c < input.nSecondStageCons(scen) && allRowsFeasible; c++)
			{
				int brokenDirection=rootSolver.isRowFeasible(c, scen, roundedSolution);
				//cout<<"Trying to fix constraint "<<c<<" s"<<scen<<" which is broken on the side "<<brokenDirection<<endl;
		    		
				if (brokenDirection!=1){
					allRowsFeasible=false;
					
					const CoinShallowPackedVector row2=rootSolver.retrieveWRow(c,scen);
			    	int nElems2=row2.getNumElements();
					const int*indices2=row2.getIndices();
					const double *elems2=row2.getElements(); 
					double bestLock=COIN_DBL_MAX;
					int varIndexOfBestLock=-1;
					int scenOfBestLock=-1;
					double objCont=COIN_DBL_MAX;//We are trying to minimize
					int directionToRound;
			
			    	for (int el=0; el<nElems2; el++){
			    		if (input.isSecondStageColInteger(scen,indices2[el]) && !isIntFeas(roundedSolution.getSecondStageVec(scen)[indices2[el]],intTol)){ //We have a possible rounding candidate
			    			if ((elems2[el]>0 && brokenDirection==-1) || (elems2[el]<0 && brokenDirection==-2) ){ ///If we have the constraint broken in the ub and the cofficient is positive, OR we have the constraint broken in the lb and the coefficient is negative we must round down
			    				if (bestLock>downLocks.getSecondStageVec(scen)[indices2[el]] ||(bestLock==downLocks.getSecondStageVec(scen)[indices2[el]] && objCont> objContribution2(roundedSolution.getSecondStageVec(scen)[indices2[el]],varObjectives.getSecondStageVec(scen)[indices2[el]],0))){
			    					bestLock=downLocks.getSecondStageVec(scen)[indices2[el]];
			    					varIndexOfBestLock=indices2[el];
			    					scenOfBestLock=scen;
			    					directionToRound=0;
			    					objCont=objContribution2(roundedSolution.getSecondStageVec(scen)[indices2[el]],varObjectives.getSecondStageVec(scen)[indices2[el]],0);
			    				}
			    			}
			    			else if ((elems2[el]<0 && brokenDirection==-1) || (elems2[el]>0 && brokenDirection==-2)){ ///If we have the constraint broken in the ub and the cofficient is negative, OR we have the constraint broken in the ub and the coefficient is positive we must round down
			    				if (bestLock>upLocks.getSecondStageVec(scen)[indices2[el]] ||(bestLock==upLocks.getSecondStageVec(scen)[indices2[el]] && objCont> objContribution2(roundedSolution.getSecondStageVec(scen)[indices2[el]],varObjectives.getSecondStageVec(scen)[indices2[el]],1))){
			    					bestLock=upLocks.getSecondStageVec(scen)[indices2[el]];
			    					varIndexOfBestLock=indices2[el];
			    					scenOfBestLock=scen;
			    					directionToRound=1;
			    					objCont=objContribution2(roundedSolution.getSecondStageVec(scen)[indices2[el]],varObjectives.getSecondStageVec(scen)[indices2[el]],1);
			    				
			    				}
			    			}
			    		}    		
			    	}
			    
					if (varIndexOfBestLock != -1 && directionToRound==0){//Round down variable
						//cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<floor(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock])<<endl;
						roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=floor(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						lb.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=floor(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						ub.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=floor(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						
						return 2;
					}
					else if (varIndexOfBestLock != -1 && directionToRound==1){//Round up variable
						//cout<<"rounded var "<<varIndexOfBestLock<<" to be "<<ceil(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock])<<endl;
						roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=ceil(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						lb.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=ceil(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						ub.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]=ceil(roundedSolution.getSecondStageVec(scenOfBestLock)[varIndexOfBestLock]);
						
						return 2;
					}

				}
			}
		}

	}

	return allRowsFeasible;
}
int findFreshSecondStageVar(denseBAVector &roundedSolution,denseBAVector & upLocks, denseBAVector &downLocks,denseBAVector &lb,denseBAVector &ub){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	 BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    const denseBAVector varObjectives = rootSolver.getVarObjective();
    
	int bestLockIndex=-1;
	int bestLockScen=-1;
	double maxLock=COIN_DBL_MIN;
	int roundDirection;
	double objCont=COIN_DBL_MAX;//We are trying to minimize
	//cout<<"FF2"<<endl;	
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int v=0; v<input.nSecondStageVars(scen); v++){
				if(input.isSecondStageColInteger(scen,v) && !isIntFeas(roundedSolution.getSecondStageVec(scen)[v],intTol)){
					//We have a possible candidate
					if(maxLock<upLocks.getSecondStageVec(scen)[v] ||(maxLock==upLocks.getSecondStageVec(scen)[v] && objCont> objContribution2(roundedSolution.getSecondStageVec(scen)[v],varObjectives.getSecondStageVec(scen)[v],0))){
						maxLock=upLocks.getSecondStageVec(scen)[v];
						bestLockIndex=v;
						roundDirection=0;
						bestLockScen=scen;
						objCont=objContribution2(roundedSolution.getSecondStageVec(scen)[v],varObjectives.getSecondStageVec(scen)[v],0);
					}
					if(maxLock<downLocks.getSecondStageVec(scen)[v]||(maxLock==downLocks.getSecondStageVec(scen)[v] && objCont> objContribution2(roundedSolution.getSecondStageVec(scen)[v],varObjectives.getSecondStageVec(scen)[v],1))){
						maxLock=downLocks.getSecondStageVec(scen)[v];
						bestLockIndex=v;
						roundDirection=1;
						bestLockScen=scen;
						objCont=objContribution2(roundedSolution.getSecondStageVec(scen)[v],varObjectives.getSecondStageVec(scen)[v],1);
					}
				}
			}
		}
	}
	if (bestLockIndex!=-1){
		if (roundDirection==1) {
			//cout<<"rounded var "<<bestLockIndex<<" s"<<bestLockScen<<" "<<maxLock<<" to be "<<ceil(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex])<<endl;	
			roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			lb.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			ub.getSecondStageVec(bestLockScen)[bestLockIndex]=ceil(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			
			return 1;
		}
		else{
			//cout<<"rounded var "<<bestLockIndex<<" s"<<bestLockScen<<" "<<maxLock<<" to be "<<floor(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex])<<endl;
			roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			lb.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			ub.getSecondStageVec(bestLockScen)[bestLockIndex]=floor(roundedSolution.getSecondStageVec(bestLockScen)[bestLockIndex]);
			
			return 1;
		} 
	}
	
	return 0;
}

bool BBSMPSHeuristicLockRounding::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB){
	
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Fix and Dive heuristic.";
	timesCalled++;
	//cout<<"ch1"<<endl;
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	
	const BADimensionsSlacks &originalDimensions= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
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
    //cout<<"ch2"<<endl;
    generateLocks2(upLocks,downLocks);
    //cout<<"ch3"<<endl;
    

    denseBAVector roundedSolution(LPRelaxationSolution);

    bool isLPFeasible = isLPIntFeas(roundedSolution);
    denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	node->getAllBranchingInformation(lb,ub);
	//cout<<"ch4"<<endl;

	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	node->reconstructWarmStartState(ps);
	rootSolver.setStates(ps);
//cout<<"ch5"<<endl;
    while (!isLPFeasible){
    	////cout<<"ITERATIOOOOOOOOOON!!! "<<endl;
	    bool done=false;
	    //cout<<"ch6"<<endl;
	   	while (!done){

	    	int result= findAndFixFirstStageConstraint(roundedSolution,upLocks,downLocks,lb,ub);
	    	//cout<<"ch7"<<endl;
	    	if (result<2){
	    		//Then we have failed at recovering or it is all good. Either case, check for more integer vars.
	    		int fresh1VarRes=findFreshFirstStageVar(roundedSolution,upLocks,downLocks,lb,ub);
	    		if (fresh1VarRes==0){
	    			done=true;
	    		}
	    	}
	    	//cout<<"ch8"<<endl;
	    }

	    
	    rootSolver.setLB(lb);
		rootSolver.setUB(ub);

		
		rootSolver.commitStates();


		//Check if feasible
		rootSolver.go();

		/* Check solver status for infeasibility/optimality */
		solverState lpStatus = rootSolver.getStatus();
		bool otherThanOptimal = (Optimal != lpStatus); 
		if (otherThanOptimal) return false;

	    roundedSolution=(rootSolver.getPrimalSolution());

	    //cout<<"ch9"<<endl;
		done=false;

	    while (!done){
	    	//cout<<"ch10"<<endl;
			//Then we have failed at finding a new first stage variable. Let's look for second stage consts for fixing
			int result2= findAndFixSecondStageConstraint(roundedSolution,upLocks,downLocks,lb,ub);
			//cout<<"ch11"<<endl;
			////cout<<"We returned "<<result2<<endl;
			if (result2<2){
				//Then we have failed at recovering or it is all good. Either case, check for more 2nd stage integer vars
				int fresh2VarRes=findFreshSecondStageVar(roundedSolution,upLocks,downLocks,lb,ub);
				if (fresh2VarRes==0)done=true;
			}
			//cout<<"ch12"<<endl;
		}

		rootSolver.setLB(lb);
		rootSolver.setUB(ub);
		rootSolver.go();

		/* Check solver status for infeasibility/optimality */
		lpStatus = rootSolver.getStatus();
		otherThanOptimal = (Optimal != lpStatus); 

		int anyOtherThanOptimal=0;
	     MPI_Allreduce(&anyOtherThanOptimal, &otherThanOptimal, 1, MPI_INT,  MPI_MAX, ctx.comm());

		if (anyOtherThanOptimal>0) return false;
		
	    roundedSolution=(rootSolver.getPrimalSolution());

	    isLPFeasible = isLPIntFeas(roundedSolution);
	    //cout<<"ch13"<<endl;
	    ////cout<<"Are we LP Feasible "<<isLPFeasible<<endl;
	}

	
	solverState lpStatus = rootSolver.getStatus();
	bool otherThanOptimal = (Optimal != lpStatus); 
	
	if(!otherThanOptimal){
		denseBAVector solVector=rootSolver.getPrimalSolution();
		solution=BBSMPSSolution(solVector,rootSolver.getObjective());
		
	}
	//cout<<"ch14"<<endl;
	//return if success
	bool success= (!otherThanOptimal && rootSolver.getObjective()<objUB);
	timesSuccessful+=success;

		
	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicLockRounding::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	if (node->getNodeDepth()<10)return true;
	int nodeDepth=node->getNodeDepth();
	assert(nodeDepth>=0);
	return (nodeDepth%15==0);

}