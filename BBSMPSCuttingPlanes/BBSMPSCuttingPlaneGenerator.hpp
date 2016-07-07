// ----------------------------------------------------------------------------
/**
   File: BBSMPSCuttingPlaneGenerator.hpp

   Description: Base virtual class of a cutting plane generator.

   Limitations: An object of this class should not be instantiated, it should always
   be subclassed. It features a virtual function called "shouldItRun", which is intended
   to be a more complex boolean function to determine if the cut should be generated.

*/
// ----------------------------------------------------------------------------

#ifndef BBSMPSCUTTINGPLANEGENERATOR_H
#define BBSMPSCUTTINGPLANEGENERATOR_H

#include "BAData.hpp"
#include <cassert> // C-style assertions
#include "BBSMPSNode.hpp"
#include <string>
#include "BBSMPSLogging.hpp"
#include "BBSMPSSolution.hpp"
#include "BBSMPSCuttingPlane.hpp"

class BBSMPSCuttingPlaneGenerator {
public:
	BBSMPSCuttingPlaneGenerator(const char *_name);
	~BBSMPSCuttingPlaneGenerator();
	virtual bool shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){return true;};
	virtual bool generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){std::cout<<"Well, this is an error...\n";};

	virtual void printStatistics();
private:
		std::string name;
protected:
	int cuttingPlanesGenerated;


};


#endif

