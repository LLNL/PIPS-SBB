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
// ----------------------------------------------------------------------------
/**
   File: BBSMPSPseudoCostBranchingRule.hpp

   Description: Branching rule that branches on the maximum fractional variable first.


*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSPSEUDOCOSTBRANCHINGRULE_H
#define BBSMPSPSEUDOCOSTBRANCHINGRULE_H

#include "BBSMPSBranchingRule.hpp"
#include "BBSMPSUtils.hpp"
#include "BBSMPSSolver.hpp"
#include "BBSMPSLogging.hpp"

class BBSMPSPseudoCostBranchingRule: public BBSMPSBranchingRule {

public:

   	virtual bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln);

   	BBSMPSPseudoCostBranchingRule(int priority);


private:

	// Auxiliary functions for branching
   	int getFirstStageMinIntInfeasCol( const denseBAVector& primalSoln,  SMPSInput& input);

   	int getFirstStageMaxFint( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getFirstStageMaxFracPartCol( const denseBAVector& primalSoln,  SMPSInput& input);

	// NOTE: MPI standard requires passing ints, not bools
   	int isFirstStageIntFeas( const denseBAVector& primalSoln,  SMPSInput& input) ;

   	void branchOnFirstStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln,  SMPSInput& input) ;

   	int getSecondStageMinIntInfeasCol( const denseBAVector& primalSoln, int scen,  SMPSInput& input);

   	int isSecondStageIntFeas( const denseBAVector& primalSoln, int scen, SMPSInput & input) ;

   	void branchOnSecondStage(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln,  SMPSInput& input,BAContext &ctx,int mype) ;

    void initializeVariable(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, double lpRelaxationObjValue, denseBAVector &lb, denseBAVector &ub, int scen, int col);

    bool performRoundOfFirstStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue);

    bool performRoundOfSecondStageInitializations(const denseBAVector &nodeRelaxation,BAFlagVector<variableState> &warmstart, denseBAVector &lb, denseBAVector &ub, double lpRelaxationObjValue);

    denseBAVector downPseudoCost;
    denseBAVector upPseudoCost;
    denseBAVector downBranchingHistory;
    denseBAVector upBranchingHistory;
    bool everythingFirstStageInitialized;
    bool everythingSecondStageInitialized;
    int reliabilityFactor;

};

#endif

