#include "BBSMPSCuttingPlane.hpp"

using namespace std;
int BBSMPSCuttingPlane::planeCounter=0;

BBSMPSCuttingPlane::BBSMPSCuttingPlane(double _lb, double _ub, denseBAVector &_expr){
	lb=_lb;
	ub=_ub;
	dExpr=denseBAVector(_expr);
	isExpressionDense=true;
	uid=(++planeCounter);
}

BBSMPSCuttingPlane::BBSMPSCuttingPlane(double _lb, double _ub, sparseBAVector &_expr){
	lb=_lb;
	ub=_ub;
	sExpr=sparseBAVector(_expr);
	isExpressionDense=false;
	uid=(++planeCounter);
}

BBSMPSCuttingPlane::BBSMPSCuttingPlane(const BBSMPSCuttingPlane& p){
	lb=p.lb;
	ub=p.ub;
	if (p.isExpressionDense) dExpr=denseBAVector(p.dExpr);
	else sExpr=sparseBAVector(p.sExpr);
	isExpressionDense=p.isExpressionDense;
}


BBSMPSCuttingPlane::BBSMPSCuttingPlane(){};

BBSMPSCuttingPlane::~BBSMPSCuttingPlane(){};


bool BBSMPSCuttingPlane::applyCrossScenarioCuttingPlane(){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();

	int originalNumVars1= dimsSlacks.inner.numFirstStageVars();
	int originalTotalSize= dimsSlacks.numFirstStageVars();

	//Create as many new vars as scenarios
	int numScenarios=input.nScenarios();
 	int firstIndex=-1;
 	if(numScenarios>0){
 		firstIndex=rootSolver.addFirstStageColumn(0,COIN_DBL_MAX,1);
 		for(int i=1; i< numScenarios; i++)rootSolver.addFirstStageColumn(0,COIN_DBL_MAX,1);
 	}


	//Create new vectors for first phase: copy the elements, set to ones the new variables
 	vector<double> Acons(dimsSlacks.inner.numFirstStageVars(),0);
 	if(isExpressionDense){
 		const denseVector &cutFirstStageVec = dExpr.getFirstStageVec();
 		for(int i=0; i< originalNumVars1; i++)Acons[i]=cutFirstStageVec[i];
 	}
 	else{
 		CoinIndexedVector &v1 = sExpr.getFirstStageVec().v;
 		int nExprElems = v1.getNumElements();
 		for (int i=0; i< nExprElems; i++){
 			Acons[v1.getIndices()[i]]=v1.denseVector()[i];
 		}
 	}
 	for(int i=originalNumVars1; i< Acons.size(); i++)Acons[i]=1;

	//At this point, the A vector is ready
 	rootSolver.addFirstStageRow(Acons,lb,ub);
	//Create new vectors for each of the scenarios
	//For each scenario:
	for (int scen=0; scen<numScenarios; scen++){
	 	if(ctx.assignedScenario(scen)){

			//Generate T with ones on the variable that is local
			vector<double> Tcons(dimsSlacks.inner.numFirstStageVars(),0);

			Tcons[firstIndex+scen]=-1;

			//Generate W with the elements belonging to the second expression
			int numVars2= dimsSlacks.inner.numSecondStageVars(scen);

			vector<double> Wcons(numVars2,0);
			if(isExpressionDense){
		 		const denseVector &cutSecondStageVec = dExpr.getSecondStageVec(scen);
				for(int i=0; i< numVars2; i++)Wcons[i]=cutSecondStageVec[i];
		 	}
		 	else{
		 		CoinIndexedVector &v2 = sExpr.getSecondStageVec(scen).v;
		 		int nExprElems2 = v2.getNumElements();
		 		for (int i=0; i< nExprElems2; i++){
		 			Wcons[v2.getIndices()[i]]=v2.denseVector()[i];
		 		}
		 	}

			//Set bounds to zero and add new constraint.
			rootSolver.addRow(Tcons,Wcons,scen,0,0);

	 		}

 	}

	return true;
}

bool BBSMPSCuttingPlane::applySingleScenarioCuttingPlane(){

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int originalNumVars1= dimsSlacks.inner.numFirstStageVars();
	int originalTotalSize= dimsSlacks.numFirstStageVars();
	int numScenarios=input.nScenarios();
	int mype=BBSMPSSolver::instance()->getMype();


	for (int scen=0; scen<numScenarios; scen++){
	 	if(ctx.assignedScenario(scen)){
	 		if(sExpr.getSecondStageVec(scen).v.getNumElements()>0){
	 			vector<double> Tcons(dimsSlacks.inner.numFirstStageVars(),0);

	 			CoinIndexedVector &v1 = sExpr.getFirstStageVec().v;
		 		int nExprElems = v1.getNumElements();
		 		for (int i=0; i< nExprElems; i++){
		 			Tcons[v1.getIndices()[i]]=v1.denseVector()[v1.getIndices()[i]];
		 		}

				//Generate W with the elements belonging to the second expression
				int numVars2= dimsSlacks.inner.numSecondStageVars(scen);
				vector<double> Wcons(numVars2,0);
		 		CoinIndexedVector &v2 = sExpr.getSecondStageVec(scen).v;
		 		int nExprElems2 = v2.getNumElements();
		 		for (int i=0; i< nExprElems2; i++){
		 			Wcons[v2.getIndices()[i]]=v2.denseVector()[v2.getIndices()[i]];
		 		}

			 	//Set bounds to zero and add new constraint.
				rootSolver.addRow(Tcons,Wcons,scen,lb,ub);
				return true;
	 		}
	 	}
	 }
	 if(sExpr.getFirstStageVec().v.getNumElements()>0){
	 	vector<double> Acons(dimsSlacks.inner.numFirstStageVars(),0);
 		CoinIndexedVector &v1 = sExpr.getFirstStageVec().v;
 		int nExprElems = v1.getNumElements();
 		for (int i=0; i< nExprElems; i++){
 			Acons[v1.getIndices()[i]]=v1.denseVector()[v1.getIndices()[i]];
 		}
 		rootSolver.addFirstStageRow(Acons,lb,ub);
		return true;
		}


	return true;
}
bool BBSMPSCuttingPlane::applyCuttingPlane(){

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int mype=BBSMPSSolver::instance()->getMype();
	if (isExpressionDense){
		applyCrossScenarioCuttingPlane();
	}
	else{
		int nScensInExpr=0;
		int numScenarios=input.nScenarios();
		for (int scen=0; scen<numScenarios; scen++){
		 	if(ctx.assignedScenario(scen)){
		 		nScensInExpr+=(sExpr.getSecondStageVec(scen).v.getNumElements()>0);

		 	}
		 }
		MPI_Allreduce(MPI_IN_PLACE,&nScensInExpr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		//nScensInExpr+= (sExpr.getFirstStageVec().v.getNumElements()>0);
		//Do crossScenario cut
		if (nScensInExpr>1){
			applyCrossScenarioCuttingPlane();
		}
		else {
			//Decide on scenario and apply only if we own it.
			applySingleScenarioCuttingPlane();
		}
	}




}


