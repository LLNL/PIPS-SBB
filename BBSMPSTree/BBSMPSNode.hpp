// ----------------------------------------------------------------------------
/**
   File: BBSMPSNode.hpp

   Description: Header file for the node class. The node essentially contains
                the information relative to a specific point in the BB tree:
                - Hotstart for the solver
                - Branching information
                - (Future) cut information

   Limitations: Since a node only contains incremental information, access to 
                the node hierarchy is required in order to recover cuts and branching
                info. Make sure to link to parents properly upon instantiation.

*/ 
// ----------------------------------------------------------------------------


#ifndef BBSMPSNODE_H
#define BBSMPSNODE_H

#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "BBSMPSBranchingInfo.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert> // C-style assertions


class BBSMPSNode {
public:

  // Construct node from {lower, upper} bounds, obj val of parent node
  BBSMPSNode(double _objectiveValue,const BAFlagVector<variableState> &states);

  BBSMPSNode(BBSMPSNode* parent,std::vector<BBSMPSBranchingInfo>& bInfos);

  // Add copy constructor so that priority_queue can use it for
  // instantiation because denseBAVector and BAFlagVector do not have
  // assignment operators or copy constructors.
  BBSMPSNode(const BBSMPSNode &sourceNode);

  ~BBSMPSNode();

  // Overload assignment operator so that priority_queue can use it
  // for instantiation because denseBAVector and BAFlagVector do not
  // have assignment operators or copy constructors. Satisfies
  // "Rule of 3". (For C++11, follow "Rule of 5".)
  BBSMPSNode& operator=(const BBSMPSNode& sourceNode) {

    /* Diagnostic code: delete when no longer needed */
    /* Also assumes communicator is MPI_COMM_WORLD. */
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    //    if (0 == mype) cout << "Calling copy assignment operator!\n";
    // Check for self-assignment
    if (this == &sourceNode) {

      //if (0 == mype) cout << "Calling self-assignment branch!\n";
      return *this;
    }

    // Copy-assign each member individually
    warmStartState = sourceNode.warmStartState;
    parent = sourceNode.parent;
    objectiveValue = sourceNode.objectiveValue;
    childrenAlive = sourceNode.childrenAlive;
    branchingInfos = sourceNode.branchingInfos;
    nodeNumber=sourceNode.nodeNumber;
    nodeDepth=sourceNode.nodeDepth;


    // Return existing object for chaining.
    //if (0 == mype) cout << "Exiting copy assignment operator!\n";
    return *this;
  }

  double getParentObjective() const;

  double getObjective() const;

  void setObjective(double lb);

  void setParentPtr(BBSMPSNode* parent);

  BBSMPSNode* getParentPtr();

  void addBranchingInformation(BBSMPSBranchingInfo& bi);

  void incrementAliveChildren();

  void decrementAliveChildren();

  void eliminate();
  
  void setWarmStartState(BAFlagVector<variableState> &state);

  void getWarmStartState(BAFlagVector<variableState> &state);

  void getAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector);

  void getAllBranchingInformation(denseBAVector &lb,denseBAVector &ub);

  int getNodeNumber() const;

  void setNodeDepth(int depth);

  int getNodeDepth();

  static void initializeNodeCounter();


private:

  //Class variable used to assign node numbers upon instantiation
  static int nodeCounter;

  int nodeNumber;

  //Depth of the node within the BB tree
  int nodeDepth;

  //Pointer to the node's parent
  BBSMPSNode *parent;

  //Child reference counter
  int childrenAlive;

  //Incremental branching information relative to this node
  std::vector<BBSMPSBranchingInfo> branchingInfos;

  // variable states for warm start information; each index is
  // one of {Basic, AtLower, AtUpper}
  BAFlagVector<variableState> warmStartState;

  // objective function of the actual node (not its parent)
  double objectiveValue; 


  // TODO: Local cut objects; Global cuts are stored in the B&B tree.

  void auxCopyAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector);

  int getBranchingInfoSize();

  BBSMPSNode(); // Disallow default constructor

};

#endif