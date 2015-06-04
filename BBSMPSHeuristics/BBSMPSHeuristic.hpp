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

class BBSMPSHeuristic {
public:
	BBSMPSHeuristic(int offset, int depth,  const char *_name);
	~BBSMPSHeuristic();
	virtual bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){return true;};
	bool checkPeriodicity(BBSMPSNode* node);
	virtual bool runHeuristic(BBSMPSNode* node, denseBAVector &LPRelaxationSolution, denseBAVector &solution, double &solObjValue){std::cout<<"Well, this is an error...\n";};

	virtual void printStatistics();
private:
		int offset;
		int depth;
		std::string name;
protected:
		int timesCalled;
		int timesSuccessful;

};


#endif

