/*########################################################################
Copyright (c) 2014-2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

Created by Geoffrey Oxberry (oxberry1@llnl.gov, goxberry@gmail.com),
Lluis-Miquel Munguia Conejero (lluis.munguia@gatech.edu), and Deepak
Rajan (rajan3@llnl.gov). LLNL-CODE-699387. All rights reserved.

This file is part of PIPS-SBB. For details, see
https://github.com/llnl/PIPS-SBB.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (as
published by the Free Software Foundation) version 2.1, February 1999.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
conditions of the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
########################################################################*/
#include "BBSMPSHeuristicCrossover.hpp"

using namespace std;


bool BBSMPSHeuristicCrossover::runHeuristic(BBSMPSNode* node, denseBAVector &nodeSolution){

	int mype=BBSMPSSolver::instance()->getMype();

	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Performing the Crossover heuristic.";

	int originalSolutionPoolSize=BBSMPSSolver::instance()->getSolPoolSize();
	double objUB=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

	bool didRun=false;
	//Steps for the heuristic
	double startTimeStamp = MPI_Wtime();

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

	int nSols=BBSMPSSolver::instance()->getSolPoolSize();

	vector< std::pair< int, int> > crossoversToDo;

	for (int i=0; i< nSols; i++){
		const BBSMPSSolution sol1=BBSMPSSolver::instance()->getSoln(i);

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

		//Run
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
		BBSMPSHeuristicLockRounding *hr= new BBSMPSHeuristicLockRounding(1,15,"LockRounding");
		bb.loadLPHeuristic(hr);

		//Add time/node limit
		bb.setNodeLimit(nodeLim);
		double bestUB=BBSMPSSolver::instance()->getSoln(0).getObjValue();

		bb.setLB(node->getObjective());
		bb.setUB(bestUB);
		//Run


		bb.branchAndBound();
		didRun=true;

		seenCrossovers[crossoversToDo[i]]=1;


	}
	double objUB2=COIN_DBL_MAX;
	if (BBSMPSSolver::instance()->getSolPoolSize()>0)objUB2=BBSMPSSolver::instance()->getSoln(0).getObjValue();


		//Retrieve best solution and return
	bool success=(objUB!=objUB2);

	timesCalled++;
	timesSuccessful+=(success);
	cumulativeTime+=(MPI_Wtime()-startTimeStamp);
	return success&didRun;

}

bool BBSMPSHeuristicCrossover::shouldItRun(BBSMPSNode* node, denseBAVector &nodeSolution){
	return true;

}
