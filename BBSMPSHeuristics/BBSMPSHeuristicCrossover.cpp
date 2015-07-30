#include "BBSMPSHeuristicCrossover.hpp"

using namespace std;


bool BBSMPSHeuristicCrossover::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	
	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	int nSols=BBSMPSSolver::instance()->getSolPoolSize();

	vector< std::pair< int, int> > crossoversToDo;
	
	for (int i=0; i< nSols; i++){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSoln(i);
		cout<<"In the queue we have solution "<<sol1.getSolNumber()<<endl;
		for (int j=i+1; j< nSols; j++){
			const BBSMPSSolution sol2=BBSMPSSolver::instance()->getSoln(j);
			std::pair<int,int> index(sol1.getSolNumber(),sol2.getSolNumber());
			if (seenCrossovers.count(index)==0){


				//Add To Crossovers list
				crossoversToDo.push_back(index);

			}
		}
	}

	for (int i=0; i< crossoversToDo.size(); i++){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSolnBySolNumber(crossoversToDo[i].first);
		const BBSMPSSolution sol2=BBSMPSSolver::instance()->getSolnBySolNumber(crossoversToDo[i].second);
			
		int count=0;
		int totalVars=0;	
		
		denseBAVector solutionVector1;
		sol1.getSolutionVector(solutionVector1);

		denseBAVector solutionVector2;
		sol2.getSolutionVector(solutionVector2);

		//Find Differences
		//Create an empty branching info vector
		vector<BBSMPSBranchingInfo> bInfos;

		//Fix variables
		for (int col = 0; col < input.nFirstStageVars(); col++)
		{	
			if(input.isFirstStageColInteger(col)){
				double sol1Val=solutionVector1.getFirstStageVec()[col];
				double sol2Val=solutionVector2.getFirstStageVec()[col];
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
						double sol2Val=solutionVector2.getSecondStageVec(scen)[col];
						if (sol1Val==sol2Val){
							bInfos.push_back(BBSMPSBranchingInfo(col,sol1Val,'E',2,scen));
							count++;
						}
						totalVars++;
					}

				}
			}
		}

		cout<<"Running a crossover between "<<sol1.getSolNumber()<< " and "<<sol2.getSolNumber()<<" fixes "<<count<<" out of "<<totalVars<<endl;
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
		bb.loadSimpleHeuristics();
		
		//Add time/node limit
		bb.setNodeLimit(nodeLim);
		double bestUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();
	
		bb.setLB(node->getObjective());
		bb.setUB(bestUB);
		//Run
		
		bb.branchAndBound();

		
		seenCrossovers[crossoversToDo[i]]=1;
			
		
	}


		//Retrieve best solution and return
	bool success=(originalSolutionPoolSize!=BBSMPSSolver::instance()->getSolPoolSize());
	timesCalled++;
	timesSuccessful+=(success);

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return false;

/*

	  const BBSMPSSolution &getSoln(int index);
  int getSolPoolSize();



		int count=0;
	int totalCount=0;
	
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
	//BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	//rootNode.setWarmStartState(ps);

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
*/
}

bool BBSMPSHeuristicCrossover::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	return true;

}