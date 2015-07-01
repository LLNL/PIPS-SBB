#include "BBSMPSHeuristicRounding.hpp"

using namespace std;


bool BBSMPSHeuristicRounding::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, BBSMPSSolution &solution, double objUB){
	
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	timesCalled++;
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the simple rounding heuristic.";
		

	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
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

	BAFlagVector<variableState> ps;
	node->getWarmStartState(ps);
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
	
	if(!otherThanOptimal){
		denseBAVector solVector=rootSolver.getPrimalSolution();
		solution=BBSMPSSolution(solVector,rootSolver.getObjective());
		
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
	return ((numberOfFractionalVariables*100/nIntVars)<25 );

}