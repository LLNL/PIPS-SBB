#include "BBSMPSHeuristicRENS.hpp"

using namespace std;


bool BBSMPSHeuristicRENS::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	
	//Steps for the heuristic


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
			double solValue=nodeSolution.getFirstStageVec()[col];
			if (isIntFeas(solValue,intTol)){//Then we fix the variable
				bInfos.push_back(BBSMPSBranchingInfo(col,solValue,'E',1));
				count++;
			}
			else {
				double lb=floor(solValue);
				double ub=ceil(solValue);
				bInfos.push_back(BBSMPSBranchingInfo(col,lb,'L',1));
				bInfos.push_back(BBSMPSBranchingInfo(col,ub,'U',1));
				
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
					double solValue=nodeSolution.getSecondStageVec(scen)[col];
					if (isIntFeas(solValue,intTol)){//Then we fix the variable
						bInfos.push_back(BBSMPSBranchingInfo(col,solValue,'E',2,scen));
						count++;
					}
					else {
						double lb=floor(solValue);
						double ub=ceil(solValue);
						bInfos.push_back(BBSMPSBranchingInfo(col,lb,'L',2,scen));
						bInfos.push_back(BBSMPSBranchingInfo(col,ub,'U',2,scen));
						
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
	//bb.setNodeLimit(2);
	//Run
	bb.setNodeLimit(nodeLim);
	bb.branchAndBound();

	//Retrieve best solution and return
	bool success=bb.retrieveBestSolution(solution);
	timesCalled++;
	timesSuccessful+=(success);


	return success;

}

bool BBSMPSHeuristicRENS::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	int numberOfFreeVars=0;
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	for (int col = 0; col < input.nFirstStageVars(); col++)
	{	
		if(input.isFirstStageColInteger(col)){
			double solValue=nodeSolution.getFirstStageVec()[col];
			if (!isIntFeas(solValue,intTol)){//Then we fix the variable
				numberOfFreeVars++;
			}
			
			
		}
	}

	return (numberOfFreeVars<20);

}