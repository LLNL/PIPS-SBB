#include "BBSMPSHeuristicRounding.hpp"

using namespace std;


// Outputs solver status:
void outputLPStatus2(solverState lpStatus) {
	cout << "PIPS-S has returned ";
	string status;
	switch(lpStatus) {
		case Uninitialized:
		status = "Uninitialized";
		break;
		case LoadedFromFile:
		status = "LoadedFromFile";
		break;
		case Initialized:
		status = "Initialized";
		break;
		case PrimalFeasible:
		status = "PrimalFeasible";
		break;
		case DualFeasible:
		status = "DualFeasible";
		break;
		case Optimal:
		status = "Optimal";
		break;
		case ProvenUnbounded:
		status = "ProvenUnbounded";
		break;
		case ProvenInfeasible:
		status = "ProvenInfeasible";
		break;
		case Stopped:
		status = "Stopped";
		break;
	}
	cout << status << endl;
}


bool BBSMPSHeuristicRounding::runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, denseBAVector &solution){
	
	int mype=BBSMPSSolver::instance()->getMype();
	
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the simple rounding heuristic.";
		

	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());
	//Apply Simple rounding to 1st stage vars 
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		lb.getFirstStageVec()[col]=roundToNearestInteger(LPRelaxationSolution.getFirstStageVec()[col]);
		ub.getFirstStageVec()[col]=roundToNearestInteger(LPRelaxationSolution.getFirstStageVec()[col]);

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
	outputLPStatus2(lpStatus);
	if (otherThanOptimal) return false;

	denseBAVector primalSoln(rootSolver.getPrimalSolution());
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	//Apply simple rounding to 2nd stage vars
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				lb.getSecondStageVec(scen)[col]=roundToNearestInteger(primalSoln.getSecondStageVec(scen)[col]);
				ub.getSecondStageVec(scen)[col]=roundToNearestInteger(primalSoln.getSecondStageVec(scen)[col]);

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
	outputLPStatus2(lpStatus);
	otherThanOptimal = (Optimal != lpStatus); 
	
	if(!otherThanOptimal){
		solution=rootSolver.getPrimalSolution();
	}
	//return if success
	timesCalled++;
	timesSuccessful+=(!otherThanOptimal);

	if (0 == mype && !otherThanOptimal) BBSMPS_ALG_LOG_SEV(info) << "The simple rounding heuristic was successful.";
		

	return !otherThanOptimal;

}

bool BBSMPSHeuristicRounding::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	
	int numberOfFractionalVariables=0;
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		numberOfFractionalVariables+=(!isIntFeas(LPRelaxationSolution.getFirstStageVec()[col],intTol));
		
	}
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				numberOfFractionalVariables+=(!isIntFeas(roundToNearestInteger(LPRelaxationSolution	.getSecondStageVec(scen)[col]),intTol));
				
			}
		}
	}

	return numberOfFractionalVariables<2;

}