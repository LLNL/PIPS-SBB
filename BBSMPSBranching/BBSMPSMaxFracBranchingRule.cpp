#include "BBSMPSMaxFracBranchingRule.hpp"

using namespace std;
int BBSMPSMaxFracBranchingRule::getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input) {
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


int BBSMPSMaxFracBranchingRule::getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input){
	int col;

	double maxFracPart = 0;
	int maxFracPartCol = -1;

    // Return index of integer variable with largest fractional part
	for (col = 0; col < input.nFirstStageVars(); col++)
	{
		bool isColInteger = input.isFirstStageColInteger(col);
		double colFracPart = fracPart(primalSoln.getFirstStageVec()[col]);

		if(isColInteger && (colFracPart > maxFracPart) ) {
			maxFracPartCol = col;
			maxFracPart = colFracPart;
		}
	}

	if (maxFracPart <= intTol) return -1;
	return maxFracPartCol;

}

	// NOTE: MPI standard requires passing ints, not bools
int BBSMPSMaxFracBranchingRule::isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) {
	return (getFirstStageMinIntInfeasCol(primalSoln,input) == -1);
}


void BBSMPSMaxFracBranchingRule::branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input) {

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

int BBSMPSMaxFracBranchingRule::getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,   SMPSInput& input) {
	
	int col;
	for (col = 0; col < input.nSecondStageVars(scen); col++)
	{
		bool isColInteger = input.isSecondStageColInteger(scen, col);
		bool isValInteger = isIntFeas(primalSoln.getSecondStageVec(scen)[col], intTol);

		// If the (col)th 2nd stage primal variable of the (scen)th	
		// scenario is integer, but has fractional value, return idx
		if (isColInteger && !isValInteger) return col;
	}

    // Otherwise, return -1;
	return -1;
}


int BBSMPSMaxFracBranchingRule::isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) {
	return (getSecondStageMinIntInfeasCol(primalSoln, scen, input) == -1);
}

void BBSMPSMaxFracBranchingRule::branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx, int mype) {


    // For now, find the minimum scenario number on each rank such
    // that one of its decision variables is integer infeasible.
    // Call that scenario number the branching candidate for each rank.
    // Then find the minimum branching candidate over all ranks. Branch on
    /// that scenario.

    // In the scenario number selected for branching, get the minimal
    // index of an integer infeasible variable.

	int myRankBranchScen(input.nScenarios() + 1);
    //if (0 == mype) cout << "myRankBranchScen = " << myRankBranchScen << endl;
	for (int scen = 0; scen < input.nScenarios(); scen++)
	{
		if(ctx.assignedScenario(scen)) {
			if(!isSecondStageIntFeas(primalSoln, scen,input)) {
				myRankBranchScen = scen;
				break;
			}
		}
	}

	int branchScen;
	int errorFlag = MPI_Allreduce(&myRankBranchScen,
		&branchScen,
		1,
		MPI_INT, 
		MPI_MIN,
		ctx.comm());
	if (0 == mype) BBSMPS_ALG_LOG_SEV(info) << "Branching on second stage scenario "
		<< branchScen << ".";
	if (ctx.assignedScenario(branchScen))BBSMPS_ALG_LOG_SEV(info) << "Processor " << mype << " will branch on second stage scenario "
	<< branchScen << ".";

    // Then, for that scenario number, get the minimal index of
    // an integer infeasible decision variable, and branch on that column
	std::vector<BBSMPSBranchingInfo> bInfosLeftKid;
	std::vector<BBSMPSBranchingInfo> bInfosRightKid;

	if(ctx.assignedScenario(branchScen)) {

		int branchCol = getSecondStageMinIntInfeasCol(primalSoln, branchScen,input);
		bInfosLeftKid.push_back( BBSMPSBranchingInfo(branchCol, ceil(primalSoln.getSecondStageVec(branchScen)[branchCol]), 'L', 2,branchScen));
		bInfosRightKid.push_back( BBSMPSBranchingInfo(branchCol, floor(primalSoln.getSecondStageVec(branchScen)[branchCol]), 'U', 2,branchScen));
	}

	
	BBSMPSNode *leftKidNode= new BBSMPSNode(node, bInfosLeftKid);
	BBSMPSNode *rightKidNode= new BBSMPSNode(node, bInfosRightKid);
	childNodes.push_back(leftKidNode);
	childNodes.push_back(rightKidNode);

}


bool BBSMPSMaxFracBranchingRule::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln){
	/* Branching */
	// Decide which stage to branch on:
	// If first stage decision variables not integer feasible,
	// branch on a first stage variable, go to start of loop
	timesCalled++;
	SMPSInput &input= BBSMPSSolver::instance()->getSMPSInput();
	BAContext &ctx= BBSMPSSolver::instance()->getBAContext();
	int mype=BBSMPSSolver::instance()->getMype();

	if(!isFirstStageIntFeas(primalSoln,input)) {
		branchOnFirstStage(node, childNodes, primalSoln,input);
		timesSuccessful++;
		return true;
	}

    // If we get to this point, we know that the first stage decision variables
    // are integer feasible, but the LP solution is not integer feasible, so
    // one of the second stage scenarios must not be integer feasible, and
    // one of the variables in one of those scenarios should be branched on.
	branchOnSecondStage(node,childNodes,primalSoln,input,ctx,mype);
	timesSuccessful+=(childNodes.size()>0);
	return (childNodes.size()>0);

}