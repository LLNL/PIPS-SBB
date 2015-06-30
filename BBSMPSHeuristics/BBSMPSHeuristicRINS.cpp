#include "BBSMPSHeuristicRINS.hpp"

using namespace std;


bool BBSMPSHeuristicRINS::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	
	//Steps for the heuristic
		double startTimeStamp = MPI_Wtime();
	//Retrieve relaxation

	const denseBAVector &LPrelaxation=BBSMPSSolver::instance()->getLPRelaxation();
	int mype=BBSMPSSolver::instance()->getMype();
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	int count=0;
	int totalCount=0;
	//Create an empty branching info vector
	vector<BBSMPSBranchingInfo> bInfos;
	//Every time we see a match between the LPrelaxation and the Solution create a new branching info
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			double LPRValue=LPrelaxation.getFirstStageVec()[col];
			double solValue=nodeSolution.getFirstStageVec()[col];
			if (fabs(LPRValue-solValue)<intTol&& isIntFeas(solValue,intTol)){//Then we fix the variable
				bInfos.push_back(BBSMPSBranchingInfo(col,solValue,'E',1));
				count++;
			}
			totalCount++;
		}
	}

	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				if(input.isSecondStageColInteger(scen,col)){
					double LPRValue=LPrelaxation.getSecondStageVec(scen)[col];
					double solValue=nodeSolution.getSecondStageVec(scen)[col];
					if (fabs(LPRValue-solValue)<intTol && isIntFeas(solValue,intTol)){//Then we fix the variable
						bInfos.push_back(BBSMPSBranchingInfo(col,solValue,'E',2,scen));
						count++;
					}
					totalCount++;
				}
			}
		}
	}
	//Create a node
	BBSMPSNode rootNode(NULL, bInfos);
	BAFlagVector<variableState> ps;
	node->getWarmStartState(ps);
	rootNode.setWarmStartState(ps);

	//Create a tree && Add node to tree
	BBSMPSTree bb(rootNode,COIN_DBL_MIN,objUB);
	bb.setVerbosity(false);
	//Add simple heuristics to tree
	bb.loadSimpleHeuristics();
	
	//Add time/node limit
	bb.setNodeLimit(nodeLim);
	//Run

	bb.branchAndBound();

	//Retrieve best solution and return
	bool success=bb.retrieveBestSolution(solution);
	timesCalled++;
	timesSuccessful+=(success);

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicRINS::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	
	int numberOfFreeVars=0;
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	const denseBAVector &LPrelaxation=BBSMPSSolver::instance()->getLPRelaxation();
	
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			double LPRValue=LPrelaxation.getFirstStageVec()[col];
			double solValue=nodeSolution.getFirstStageVec()[col];
			if (fabs(LPRValue-solValue)>intTol && isIntFeas(solValue,intTol)){//Then we fix the variable
				numberOfFreeVars++;
			}
			
		}
	}
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{

				if(input.isSecondStageColInteger(scen,col)){
					double LPRValue=LPrelaxation.getFirstStageVec()[col];
					double solValue=nodeSolution.getFirstStageVec()[col];
					if (fabs(LPRValue-solValue)>intTol && isIntFeas(solValue,intTol)){//Then we fix the variable
						numberOfFreeVars++;
					}
				}
			}
		}
	}

	int minCount;
	int errorFlag = MPI_Allreduce(&numberOfFreeVars,
		&minCount,
		1,
		MPI_INT, 
		MPI_MIN,
		ctx.comm());

	return (minCount<40);
}