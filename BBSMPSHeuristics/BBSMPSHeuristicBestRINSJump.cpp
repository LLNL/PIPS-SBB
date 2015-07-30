#include "BBSMPSHeuristicBestRINSJump.hpp"

using namespace std;


bool BBSMPSHeuristicBestRINSJump::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	
	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	int nSols=BBSMPSSolver::instance()->getSolPoolSize();

	vector< int > crossoversToDo;
	const denseBAVector &LPrelaxation=BBSMPSSolver::instance()->getLPRelaxation();
	if (nSols==0){
		timesCalled++;
		cumulativeTime+=(MPI_Wtime()-startTimeStamp);
		return false;
	}

	bool exit=false;
	while(!exit){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSoln(0);
		
		if (seenCrossovers.count(sol1.getSolNumber())==0){
			seenCrossovers[sol1.getSolNumber()]=1;
			//Perform run
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
					double sol2Val=LPrelaxation.getFirstStageVec()[col];
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
							double sol2Val=LPrelaxation.getSecondStageVec(scen)[col];
							if (sol1Val==sol2Val){
								bInfos.push_back(BBSMPSBranchingInfo(col,sol1Val,'E',2,scen));
								count++;
							}
							totalVars++;
						}

					}
				}
			}

			cout<<"Running a crossover between "<<sol1.getSolNumber()<<" and relaxation. fixes "<<count<<" out of "<<totalVars<<endl;
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
			BBSMPSHeuristicLockRounding *hr= new BBSMPSHeuristicLockRounding(1,1,"LockRounding");
	  		bb.loadHeuristic(hr);
			//Add time/node limit
			bb.setNodeLimit(nodeLim);
			bb.setSolLimit(1);
			double bestUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();
		
			bb.setLB(node->getObjective());
			bb.setUB(bestUB);
			//Run
			
			bb.branchAndBound();


		}
		else exit=true;


	}


		//Retrieve best solution and return
	bool success=(originalSolutionPoolSize!=BBSMPSSolver::instance()->getSolPoolSize());
	timesCalled++;
	timesSuccessful+=(success);

	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return false;

}

bool BBSMPSHeuristicBestRINSJump::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	return true;

}