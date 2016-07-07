#include "BBSMPSHeuristicRENS.hpp"

using namespace std;


bool BBSMPSHeuristicRENS::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution){
	int mype=BBSMPSSolver::instance()->getMype();

	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the RENS heuristic.";
	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	double objUB=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();


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
	BBSMPSNode* rootNode= new BBSMPSNode(NULL, bInfos);
	//BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	//rootNode.setWarmStartState(ps);
	//rootNode->copyCuttingPlanes(BBSMPSTree::getRootNode());
	//Create a tree && Add node to tree
	BBSMPSTree bb(rootNode,COIN_DBL_MIN,objUB);
	bb.setVerbosity(false);
	//Add simple heuristics to tree
	//bb.loadSimpleHeuristics();
	BBSMPSHeuristicLockRounding *hr= new BBSMPSHeuristicLockRounding(0,15,"LockRounding");
	bb.loadLPHeuristic(hr);

	//Add time/node limit
	bb.setNodeLimit(nodeLim);
	//Run

	bb.branchAndBound();
	double objUB2=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB2=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	//Retrieve best solution and return
	bool success=(objUB!=objUB2);
	timesCalled++;
	timesSuccessful+=(success);

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success;

}

bool BBSMPSHeuristicRENS::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	int numberOfFreeVars=0;
	int numberOfIntegerVars=0;
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	const denseBAVector &lb=BBSMPSSolver::instance()->getOriginalLB();
	const denseBAVector &ub=BBSMPSSolver::instance()->getOriginalUB();

	for (int col = 0; col < input.nFirstStageVars(); col++)
	{
		if(input.isFirstStageColInteger(col)){
			double solValue=nodeSolution.getFirstStageVec()[col];
			if (!isIntFeas(solValue,intTol)){//Then we free the variable
				numberOfFreeVars++;
			}
			if(ub.getFirstStageVec()[col]-lb.getFirstStageVec()[col]>1)numberOfIntegerVars++;


		}
	}

		int numberOfFreeVars2=0;
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			for (int col = 0; col < input.nSecondStageVars(scen); col++)
			{
				if(input.isSecondStageColInteger(scen,col)){
					double solValue=nodeSolution.getSecondStageVec(scen)[col];
					if (!isIntFeas(solValue,intTol)){//Then we free the variable
						numberOfFreeVars2++;

					}
					if(ub.getSecondStageVec(scen)[col]-lb.getSecondStageVec(scen)[col]>1)numberOfIntegerVars++;
				}
			}
		}
	}
	int count;
	int errorFlag = MPI_Allreduce(&numberOfFreeVars2,
		&count,
		1,
		MPI_INT,
		MPI_SUM,
		ctx.comm());
	count+=numberOfFreeVars;
	return (count<nodeLim);

}