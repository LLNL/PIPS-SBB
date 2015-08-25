// ----------------------------------------------------------------------------
/**
   File: BBSMPSBranchingRuleManager.hpp

   Description: Branching rule manager responsible for organizing branching rules by priority
   				and calling them recursively upon node branching.

   Limitations: Branching rules are only ordered by priority.

*/ 
// ----------------------------------------------------------------------------

#ifndef BBSMPSBRANCHINGRULEMANAGER_H
#define BBSMPSBRANCHINGRULEMANAGER_H

#include <set>
#include "BBSMPSBranchingRule.hpp"

struct BRPriorityOrdering
{
	bool operator() (const BBSMPSBranchingRule* lhs, const BBSMPSBranchingRule* rhs) const
	{	
		return (lhs->getPriority() > rhs->getPriority());
	}
};


class BBSMPSBranchingRuleManager {

public:
	
	BBSMPSBranchingRuleManager();
	
	~BBSMPSBranchingRuleManager();

	void addBranchingRule(BBSMPSBranchingRule *rule);

	bool branch(BBSMPSNode * node, std::vector<BBSMPSNode*> &childNodes,  const denseBAVector& primalSoln);
	
	void printStatistics();

private:

	std::multiset<BBSMPSBranchingRule*, BRPriorityOrdering> branchingRuleList;

};

#endif