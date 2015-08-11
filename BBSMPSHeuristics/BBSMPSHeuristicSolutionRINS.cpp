#include "BBSMPSHeuristicSolutionRINS.hpp"

using namespace std;


bool BBSMPSHeuristicSolutionRINS::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	
	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSoln(0);
			
	int count=0;
	int totalVars=0;	
	
	denseBAVector solutionVector1;
	sol1.getSolutionVector(solutionVector1);

	//Find Differences
	//Create an empty branching info vector
	vector<BBSMPSBranchingInfo> bInfos;

	//Fix variables
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			double sol1Val=solutionVector1.getFirstStageVec()[col];
			double sol2Val=nodeSolution.getFirstStageVec()[col];
			if (sol1Val==sol2Val){
				bInfos.push_back(BBSMPSBranchingInfo(col,sol1Val,'E',1));
				count++;
			}
			totalVars++;
		}
	}
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				if(input.isSecondStageColInteger(scen,col)){
					double sol1Val=solutionVector1.getSecondStageVec(scen)[col];
					double sol2Val=nodeSolution.getSecondStageVec(scen)[col];
					if (sol1Val==sol2Val){
						bInfos.push_back(BBSMPSBranchingInfo(col,sol1Val,'E',2,scen));
						count++;
					}
					totalVars++;
				}

			}
		}
	}	

	//Run 
	//Create a node
	BBSMPSNode rootNode(NULL, bInfos);
	//BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	//rootNode.setWarmStartState(ps);

	//Create a tree && Add node to tree
	BBSMPSTree bb(rootNode,COIN_DBL_MIN,objUB);
	bb.setVerbosity(false);
	//Add simple heuristics to tree
	//bb.loadSimpleHeuristics();
	BBSMPSHeuristicLockRounding *hr= new BBSMPSHeuristicLockRounding(1,15,"LockRounding");
	bb.loadHeuristic(hr);
	//Add time/node limit
	bb.setNodeLimit(nodeLim);
	double bestUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	bb.setLB(node->getObjective());
	bb.setUB(bestUB);
	//Run
	
	bb.branchAndBound();

	//Retrieve best solution and return
	bool success=(originalSolutionPoolSize!=BBSMPSSolver::instance()->getSolPoolSize());
	timesCalled++;
	timesSuccessful+=(success);

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return false;

}

bool BBSMPSHeuristicSolutionRINS::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	return true;

}