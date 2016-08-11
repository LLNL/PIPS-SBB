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
#include "BBSMPSPseudoCostBranchingRule.hpp"

using namespace std;

struct doubleint { double d; int i;};
double scoreFunction(double qminus, double qplus, double mu){
  double scoreval=(1-mu)*min(qminus,qplus)+(mu)*max(qminus,qplus);
	return scoreval;
}
void BBSMPSPseudoCostBranchingRule::initializeVariable(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, double lpRelaxationObjValue, denseBAVector &lb, denseBAVector &ub, int scen, int col){

	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
    int mype=BBSMPSSolver::instance()->getMype();

    //rootSolver.setStates(warmstart);

    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();



    double diffInVal;
   	int savedLB;
  	int savedUB;
  	if(ctx.assignedScenario(scen)) {

     	savedLB= lb.getVec(scen)[col];
     	savedUB= ub.getVec(scen)[col];

     	lb.getVec(scen)[col]=ceil(nodeRelaxation.getVec(scen)[col]);
     	diffInVal=lb.getVec(scen)[col]-nodeRelaxation.getVec(scen)[col];


    }
    rootSolver.setLB(lb);
    rootSolver.setUB(ub);

    rootSolver.commitStates();

    rootSolver.go();
    if(ctx.assignedScenario(scen)) {

       double lpObj = rootSolver.getObjective();
       double diffInObjective=fabs(lpObj-lpRelaxationObjValue);
       upPseudoCost.getVec(scen)[col]+=(diffInObjective/diffInVal);
       upBranchingHistory.getVec(scen)[col]++;

       lb.getVec(scen)[col]=savedLB;
       ub.getVec(scen)[col]=floor(nodeRelaxation.getVec(scen)[col]);
       diffInVal=nodeRelaxation.getVec(scen)[col]-ub.getVec(scen)[col];


    }

  	rootSolver.setLB(lb);
    rootSolver.setUB(ub);
    //rootSolver.setStates(warmstart);
    rootSolver.commitStates();

	    rootSolver.go();

    if(ctx.assignedScenario(scen)) {

       double lpObj = rootSolver.getObjective();
       double diffInObjective=fabs(lpObj-lpRelaxationObjValue);
       downPseudoCost.getVec(scen)[col]+=diffInObjective/diffInVal;
       downBranchingHistory.getVec(scen)[col]++;
       lb.getVec(scen)[col]=savedLB;
       ub.getVec(scen)[col]=savedUB;


    }

}

bool BBSMPSPseudoCostBranchingRule::performRoundOfFirstStageInitializations(const denseBAVector &sol,BAFlagVector<variableState> &warmstart,denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue){

	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    int mype=BBSMPSSolver::instance()->getMype();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();

    rootSolver.setStates(warmstart);
    bool foundAtLeastOne=false;

	  //Fix variables
     for (int col = 0; col < input.nFirstStageVars(); col++)
     {

        if(input.isFirstStageColInteger(col) && ! isIntFeas(sol.getVec(-1)[col],intTol) && (downBranchingHistory.getVec(-1)[col]<reliabilityFactor || upBranchingHistory.getVec(-1)[col]<reliabilityFactor)){

           initializeVariable(sol,warmstart, lpRelaxationObjValue, lb, ub, -1, col);

        }
        foundAtLeastOne=foundAtLeastOne || (downBranchingHistory.getVec(-1)[col]<reliabilityFactor || upBranchingHistory.getVec(-1)[col]<reliabilityFactor);

     }

    return foundAtLeastOne;


}

bool BBSMPSPseudoCostBranchingRule::performRoundOfSecondStageInitializations(const denseBAVector &sol,BAFlagVector<variableState> &warmstart,denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue){

	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
    int mype=BBSMPSSolver::instance()->getMype();
    PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();

    rootSolver.setStates(warmstart);
    bool foundAtLeastOne=false;
    bool atLeastOneVarWasFrac=false;
	//Fix variables

    for (int scen = 0; scen < input.nScenarios(); scen++)
    {
     	int owner=-1;
     	if (ctx.assignedScenario(scen))owner=mype;
     	int sharedOwner;
     	int errorFlag = MPI_Allreduce(&owner, &sharedOwner, 1, MPI_INT,  MPI_MAX, ctx.comm());

       	if (owner!=-1){
       		for (int col = 0; col < input.nSecondStageVars(scen); col++)
        	{
           		int shouldItRun=0;
          		if(input.isSecondStageColInteger(scen,col) && !isIntFeas(sol.getVec(scen)[col],intTol) && (downBranchingHistory.getVec(scen)[col]<reliabilityFactor || upBranchingHistory.getVec(scen)[col]<reliabilityFactor))shouldItRun=1;
          		MPI_Bcast( &shouldItRun, 1, MPI_INT, sharedOwner, ctx.comm());
           		if(shouldItRun){
           			initializeVariable(sol,warmstart, lpRelaxationObjValue, lb, ub, scen, col);

           		}
           		 foundAtLeastOne=foundAtLeastOne || (downBranchingHistory.getVec(scen)[col]<reliabilityFactor || upBranchingHistory.getVec(scen)[col]<reliabilityFactor);
           	}
        }
        else{
        	for (int col = 0; col < input.nSecondStageVars(scen); col++)
        	{
            	int shouldItRun=0;
            	MPI_Bcast( &shouldItRun, 1, MPI_INT, sharedOwner, ctx.comm());
            	if(shouldItRun){
           			initializeVariable(sol,warmstart, lpRelaxationObjValue, lb, ub, scen, col);
           			foundAtLeastOne=1;
           		}
           	}
        }


    }

    int foundInteger=(foundAtLeastOne);
    int globalValue;
    MPI_Allreduce(&foundInteger, &globalValue, 1, MPI_INT,  MPI_MAX, ctx.comm());
    return globalValue>0;
}


BBSMPSPseudoCostBranchingRule::BBSMPSPseudoCostBranchingRule(int priority): BBSMPSBranchingRule(priority){
	name="Max Fractional Branching Rule";
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    BAFlagVector<variableState> warmstart(BBSMPSSolver::instance()->getOriginalWarmStart());

    downPseudoCost.allocate(dimsSlacks, ctx, PrimalVector);
    upPseudoCost.allocate(dimsSlacks, ctx, PrimalVector);
    downBranchingHistory.allocate(dimsSlacks, ctx, PrimalVector);
    upBranchingHistory.allocate(dimsSlacks, ctx, PrimalVector);
    downPseudoCost.clear();
    upPseudoCost.clear();
    downBranchingHistory.clear();
    upBranchingHistory.clear();
    reliabilityFactor=4;
    everythingFirstStageInitialized=false;
    everythingSecondStageInitialized=false;
    denseBAVector sol(BBSMPSSolver::instance()->getLPRelaxation());
    denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
    denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());


    double lpRelaxationObjValue=BBSMPSSolver::instance()->getLPRelaxationObjectiveValue();

   	performRoundOfFirstStageInitializations(sol,warmstart,lb,ub, lpRelaxationObjValue);



  };

int BBSMPSPseudoCostBranchingRule::getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input) {
int col;

// Return first index of integer variable with fractional value
for (col = 0; col < input.nFirstStageVars(); col++)
{

	bool isColInteger = input.isFirstStageColInteger(col);
	bool isValInteger = isIntFeas(primalSoln.getFirstStageVec()[col], intTol); //TODO:: Tolerance is hardcoded for now

	// If the (col)th 1st stage primal variable is integer,
	// but has a fractional value, return idx
	if(isColInteger && !isValInteger) return col;
}

// Otherwise, 1st stage is integer feasible: return -1;
return -1;
}


int BBSMPSPseudoCostBranchingRule::getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input){
int col;

double maxScore = -1;
int maxCol = -1;


// Return index of integer variable with largest fractional part
for (col = 0; col < input.nFirstStageVars(); col++)
{
	bool isColInteger = input.isFirstStageColInteger(col);
	bool isValInteger = isIntFeas(primalSoln.getFirstStageVec()[col], intTol);


	if(isColInteger && !isValInteger){
		double downDiff=primalSoln.getFirstStageVec()[col]-floor(primalSoln.getFirstStageVec()[col]);
		double upDiff=ceil(primalSoln.getFirstStageVec()[col])-primalSoln.getFirstStageVec()[col];
		double downTimes=downBranchingHistory.getFirstStageVec()[col];
		double upTimes=downBranchingHistory.getFirstStageVec()[col];

		double score=scoreFunction((downPseudoCost.getFirstStageVec()[col]/downTimes)*downDiff, (upPseudoCost.getFirstStageVec()[col]/upTimes)*upDiff, 0.1667);

    if (score > maxScore)  {
			maxCol = col;
			maxScore = score;
		}
	}
}
  double downTimes=downBranchingHistory.getFirstStageVec()[maxCol];
	double upTimes=downBranchingHistory.getFirstStageVec()[maxCol];
	double downDiff=primalSoln.getFirstStageVec()[maxCol]-floor(primalSoln.getFirstStageVec()[maxCol]);
	double upDiff=ceil(primalSoln.getFirstStageVec()[maxCol])-primalSoln.getFirstStageVec()[maxCol];

//if (maxScore <= intTol) return -1;
return maxCol;

}

// NOTE: MPI standard requires passing ints, not bools
int BBSMPSPseudoCostBranchingRule::isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) {
return (getFirstStageMinIntInfeasCol(primalSoln,input) == -1);
}


void BBSMPSPseudoCostBranchingRule::branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input) {

int mype=BBSMPSSolver::instance()->getMype();
/* Branching Rule */
// For now, get minimal index of an integer infeasible variable
//int branchCol = getFirstStageMinIntInfeasCol(primalSoln);

// Get index of maximum fractional part.
int branchCol = getFirstStageMaxFracPartCol(primalSoln,input);
assert(branchCol > -1); // Should always be true if not integer feasible

if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Branching on first stage variable "<< branchCol <<".";

//Create both branching infos
std::vector<BBSMPSBranchingInfo> bInfosLeftKid;
std::vector<BBSMPSBranchingInfo> bInfosRightKid;
bInfosLeftKid.push_back( BBSMPSBranchingInfo(branchCol, ceil(primalSoln.getFirstStageVec()[branchCol]), 'L', 1));
bInfosRightKid.push_back( BBSMPSBranchingInfo(branchCol, floor(primalSoln.getFirstStageVec()[branchCol]), 'U', 1));
//Create both children

BBSMPSNode *leftKidNode= new BBSMPSNode(node, bInfosLeftKid);
BBSMPSNode *rightKidNode= new BBSMPSNode(node, bInfosRightKid);
childNodes.push_back(leftKidNode);
childNodes.push_back(rightKidNode);

}



int BBSMPSPseudoCostBranchingRule::getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,   SMPSInput& input) {


int col;
double maxScore = 0;
int maxCol = -1;
for (col = 0; col < input.nSecondStageVars(scen); col++)
{
	bool isColInteger = input.isSecondStageColInteger(scen, col);

	bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col], intTol);
	// If the (col)th 2nd stage primal variable of the (scen)th
	// scenario is integer, but has fractional value, return idx
	if (isColInteger && !isValInteger) {
		double downDiff=primalSoln.getSecondStageVec(scen)[col]-floor(primalSoln.getSecondStageVec(scen)[col]);
		double upDiff=ceil(primalSoln.getSecondStageVec(scen)[col])-primalSoln.getSecondStageVec(scen)[col];
		double downTimes=downBranchingHistory.getSecondStageVec(scen)[col];
		double upTimes=downBranchingHistory.getSecondStageVec(scen)[col];


		double score=scoreFunction((downPseudoCost.getSecondStageVec(scen)[col]/downTimes)*downDiff, (upPseudoCost.getSecondStageVec(scen)[col]/upTimes)*upDiff, 0.1667);

	 	if (score > maxScore)  {
			maxCol = col;
			maxScore = score;
		}

	}

}

//if (maxScore <= intTol) return -1;
return maxCol;
}


int BBSMPSPseudoCostBranchingRule::isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) {
return (getSecondStageMinIntInfeasCol(primalSoln, scen, input) == -1);
}

void BBSMPSPseudoCostBranchingRule::branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx, int mype) {


// For now, find the minimum scenario number on each rank such
// that one of its decision variables is integer infeasible.
// Call that scenario number the branching candidate for each rank.
// Then find the minimum branching candidate over all ranks. Branch on
/// that scenario.

// In the scenario number selected for branching, get the minimal
// index of an integer infeasible variable.

int myRankBranchScen(input.nScenarios() + 1);
double maxScore=-1;
double maxScen=-1;
double maxCol=-1;
for (int scen = 0; scen < input.nScenarios(); scen++)
{
	if(ctx.assignedScenario(scen)) {
		int col = getSecondStageMinIntInfeasCol(primalSoln, scen, input);
		if (col!= -1){
			double downDiff=primalSoln.getSecondStageVec(scen)[col]-floor(primalSoln.getSecondStageVec(scen)[col]);
			double upDiff=ceil(primalSoln.getSecondStageVec(scen)[col])-primalSoln.getSecondStageVec(scen)[col];
			double downTimes=downBranchingHistory.getSecondStageVec(scen)[col];
			double upTimes=downBranchingHistory.getSecondStageVec(scen)[col];


			double score=scoreFunction((downPseudoCost.getSecondStageVec(scen)[col]/downTimes)*downDiff, (upPseudoCost.getSecondStageVec(scen)[col]/upTimes)*upDiff, 0.1667);

			if (score>maxScore){
				maxScore=score;
				maxScen=scen;
				maxCol=col;
			}
		}
	}
}

doubleint my = { maxScore, mype }, best;
MPI_Allreduce(&my,&best,1,MPI_DOUBLE_INT,MPI_MAXLOC,ctx.comm());



// Then, for that scenario number, get the minimal index of
// an integer infeasible decision variable, and branch on that column
std::vector<BBSMPSBranchingInfo> bInfosLeftKid;
std::vector<BBSMPSBranchingInfo> bInfosRightKid;

if(best.i==mype) {

	bInfosLeftKid.push_back( BBSMPSBranchingInfo(maxCol, ceil(primalSoln.getSecondStageVec(maxScen)[maxCol]), 'L', 2,maxScen));
	bInfosRightKid.push_back( BBSMPSBranchingInfo(maxCol, floor(primalSoln.getSecondStageVec(maxScen)[maxCol]), 'U', 2,maxScen));
}


BBSMPSNode *leftKidNode= new BBSMPSNode(node, bInfosLeftKid);
BBSMPSNode *rightKidNode= new BBSMPSNode(node, bInfosRightKid);
childNodes.push_back(leftKidNode);
childNodes.push_back(rightKidNode);

}


bool BBSMPSPseudoCostBranchingRule::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln){
/* Branching */
// Decide which stage to branch on:
// If first stage decision variables not integer feasible,
// branch on a first stage variable, go to start of loop
timesCalled++;

SMPSInput &input= BBSMPSSolver::instance()->getSMPSInput();
BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
int mype=BBSMPSSolver::instance()->getMype();

if (!everythingFirstStageInitialized){

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());

	node->getAllBranchingInformation(lb,ub);
	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	double lpRelaxationObjValue=node->getObjective();
	everythingFirstStageInitialized=!performRoundOfFirstStageInitializations(primalSoln,ps,lb,ub, lpRelaxationObjValue);

}

if(!isFirstStageIntFeas(primalSoln,input)) {
	branchOnFirstStage(node, childNodes, primalSoln,input);
	timesSuccessful++;

	return true;
}



// If we get to this point, we know that the first stage decision variables
// are integer feasible, but the LP solution is not integer feasible, so
// one of the second stage scenarios must not be integer feasible, and
// one of the variables in one of those scenarios should be branched on.

if (!everythingSecondStageInitialized){

	denseBAVector lb(BBSMPSSolver::instance()->getOriginalLB());
	denseBAVector ub(BBSMPSSolver::instance()->getOriginalUB());

	node->getAllBranchingInformation(lb,ub);
	BAFlagVector<variableState> ps(BBSMPSSolver::instance()->getOriginalWarmStart());
	//node->reconstructWarmStartState(ps);
	double lpRelaxationObjValue=node->getObjective();
	everythingSecondStageInitialized=!performRoundOfSecondStageInitializations(primalSoln,ps,lb,ub, lpRelaxationObjValue);

}

branchOnSecondStage(node,childNodes,primalSoln,input,ctx,mype);
timesSuccessful+=(childNodes.size()>0);
return (childNodes.size()>0);

}
