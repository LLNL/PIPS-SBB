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
#include "BBSMPSBranchingRuleManager.hpp"

BBSMPSBranchingRuleManager::BBSMPSBranchingRuleManager(){};
BBSMPSBranchingRuleManager::~BBSMPSBranchingRuleManager(){};


void BBSMPSBranchingRuleManager::addBranchingRule(BBSMPSBranchingRule *rule){
	assert(rule!=NULL);
	branchingRuleList.insert(rule);
}

bool BBSMPSBranchingRuleManager::branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes, const denseBAVector& primalSoln){
	bool success=false;
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end() && !success; ++it){
		BBSMPSBranchingRule *br=(*it);
		success= success || br->branch(node, childNodes,  primalSoln);

	}
	return success;
}


void BBSMPSBranchingRuleManager::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"**************BRANCHING STATISTICS****************";
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end(); ++it){
		BBSMPSBranchingRule *br=(*it);
		br->printStatistics();
	}
	BBSMPS_ALG_LOG_SEV(warning)<<"**************************************************";
}

void BBSMPSBranchingRuleManager::freeResources(){
	std::multiset<BBSMPSBranchingRule*>::iterator it;
	for (it=branchingRuleList.begin(); it!=branchingRuleList.end(); ++it){
		BBSMPSBranchingRule *br=(*it);
		delete br;
	}
}
