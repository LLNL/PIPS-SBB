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
   File: BBSMPSHeuristic.hpp

   Description: Base virtual class of a primal solution heuristic.

   Limitations: An object of this class should not be instantiated, it should always
   be subclassed. A heuristic must be initialized with an offset and depth frequency.
   Additionally, it features a virtual function called "shouldItRun", which is intended
   to be a more complex boolean function to determine if the heuristic should run.

*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSHEURISTIC_H
#define BBSMPSHEURISTIC_H

#include "BAData.hpp"
#include <cassert> // C-style assertions
#include "BBSMPSNode.hpp"
#include <string>
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolution.hpp"


class BBSMPSHeuristic {
public:
	BBSMPSHeuristic(int offset, int depth,  const char *_name);
	~BBSMPSHeuristic();
	virtual bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){return true;};
	bool checkPeriodicity(BBSMPSNode* node);
	virtual bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){std::cout<<"Well, this is an error...\n";};
	double getCumulativeTime(){return cumulativeTime;};
	virtual void printStatistics();
private:

		std::string name;

protected:
		int timesCalled;
		int timesSuccessful;
		double cumulativeTime;
		int offset;
		int depth;

};


#endif

