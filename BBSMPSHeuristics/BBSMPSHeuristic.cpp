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
#include "BBSMPSHeuristic.hpp"

using namespace std;

BBSMPSHeuristic::BBSMPSHeuristic(int _offset, int _depth, const char *_name){
	offset=_offset;
	depth=_depth;
	timesCalled=0;
	name=_name;
	timesSuccessful=0;
	cumulativeTime=0;
}
BBSMPSHeuristic::~BBSMPSHeuristic(){};

bool BBSMPSHeuristic::checkPeriodicity(BBSMPSNode* node){

	int nodeDepth=node->getNodeDepth();
	assert(nodeDepth>=0);
	if (offset>nodeDepth) return false;
	nodeDepth-=offset;
	return (nodeDepth%depth==0);
}


void BBSMPSHeuristic::printStatistics(){
	BBSMPS_ALG_LOG_SEV(warning)<<"Heuristic:"<<name<<":Times Called:"<<timesCalled<<":Times successful:"<<timesSuccessful<<":Execution Time:"<<cumulativeTime;
}
