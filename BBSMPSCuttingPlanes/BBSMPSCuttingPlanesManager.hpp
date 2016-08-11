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
   File: BBSMPSCuttingPlanesManager.hpp

   Description: Cutting planes manager responsible for organizing all the cutting planes
   				and calling them recursively upon application.

   Limitations: Cutting planes are only applied when their frequency allows it and ShouldItRun()
   				returns positive.

*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSCUTTINGPLANESMANAGER_H
#define BBSMPSCUTTINGPLANESMANAGER_H

#include <set>
#include <vector>
#include "BBSMPSCuttingPlaneGenerator.hpp"
#include "BBSMPSSolver.hpp"

class BBSMPSCuttingPlanesManager {

public:
	BBSMPSCuttingPlanesManager(int maxRounds);

	~BBSMPSCuttingPlanesManager();

	void addCuttingPlaneGenerator(BBSMPSCuttingPlaneGenerator *cuttingPlaneGen);

	bool generateCuttingPlanes(BBSMPSNode* n,denseBAVector &LPRelaxationSolution);

	void printStatistics();

	void freeResources();

private:
	std::vector<BBSMPSCuttingPlaneGenerator*> cuttingPlaneGeneratorsList;
	int maxRounds;
};

#endif
