#include "BBSMPSHeuristicSolutionRINS.hpp"

using namespace std;


bool BBSMPSHeuristicSolutionRINS::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution, BBSMPSSolution &solution, double objUB){
	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	
	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();
	int mype=BBSMPSSolver::instance()->getMype();
	
	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	int nSols=BBSMPSSolver::instance()->getSolPoolSize();

	vector< int > crossoversToDo;
	const denseBAVector &LPrelaxation=BBSMPSSolver::instance()->getLPRelaxation();
	

	for (int i=0; i< nSols; i++){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSoln(i);
		cout<<"In the queue we have solution "<<sol1.getSolNumber()<<endl;
		
		if (seenCrossovers.count(sol1.getSolNumber())==0){


			//Add To Crossovers list
			crossoversToDo.push_back(sol1.getSolNumber());

			
		}
	}

	for (int i=0; i< crossoversToDo.size(); i++){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSolnBySolNumber(crossoversToDo[i]);
			
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
		BBSMPSHeuristicLockRounding *hr= new BBSMPSHeuristicLockRounding(1,15,"LockRounding");
  		bb.loadHeuristic(hr);
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

}

bool BBSMPSHeuristicSolutionRINS::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	return true;

}