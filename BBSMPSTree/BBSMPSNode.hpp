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
#include "BBSMPSCuttingPlane.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert> // C-style assertions
#include "BBSMPSLogging.hpp"

class BBSMPSNode {
public:

  // Construct node from {lower, upper} bounds, obj val of parent node
  BBSMPSNode(double _objectiveValue,const BAFlagVector<variableState> &states);
  BBSMPSNode(double _objectiveValue, const std::vector< std::pair < BAIndex, variableState > > &states);
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

  void addCuttingPlane(BBSMPSCuttingPlane *cp);

  void incrementAliveChildren();

  void decrementAliveChildren();

  void eliminate();

  void deallocateWarmStartState();

  void setIncrementalWarmStartState(std::vector< std::pair < BAIndex, variableState > > &changes);

  void reconstructWarmStartState(BAFlagVector<variableState> &state);

  void getAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector);

  void getAllBranchingInformation(denseBAVector &lb,denseBAVector &ub);

  void getAllCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector);

  void getParentNodeCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector);

  void  getGrandParentCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector);

  void getCurrentNodeCuttingPlanes(std::vector<BBSMPSCuttingPlane*> &cpVector);

  void copyCuttingPlanes(BBSMPSNode *node);

void getAllCuttingUids(std::vector<int> &uidVector);

  void getCurrentNodeCuttingPlaneUids(std::vector<int> &uidVector);
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

  std::vector< std::pair < BAIndex, variableState > > partialStartState;

  // objective function of the actual node (not its parent)
  double objectiveValue;


  //Incremental cutting plane information relative to this node
  std::vector<BBSMPSCuttingPlane*> cuttingPlanes;
  vector<int> cuttingPlaneUids;
  // TODO: Local cut objects; Global cuts are stored in the B&B tree.

  void auxCopyAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector);

  void auxCopyAllCuttingPlaneInformation(std::vector<BBSMPSCuttingPlane*> &cpVector);

  int getBranchingInfoSize();

  int getNumberOfCuttingPlanes();


  BBSMPSNode(); // Disallow default constructor

};

#endif
