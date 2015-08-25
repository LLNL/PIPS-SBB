#include "BBSMPSNode.hpp"

using namespace std;

int BBSMPSNode::nodeCounter=0;


// Construct node from {lower, upper} bounds, obj val of parent node
BBSMPSNode::BBSMPSNode(double _objectiveValue, const std::vector< std::pair < BAIndex, variableState > > &states):
partialStartState(states) {
  parent =NULL;
  objectiveValue=_objectiveValue;
  childrenAlive=0;
  nodeNumber=(++nodeCounter);
  nodeDepth=-1;
  if (nodeCounter==1) nodeDepth=0;
  
}


  // Add copy constructor so that priority_queue can use it for
  // instantiation because denseBAVector and BAFlagVector do not have
  // assignment operators or copy constructors.
BBSMPSNode::BBSMPSNode(const BBSMPSNode &sourceNode):
partialStartState(sourceNode.partialStartState),
parent(sourceNode.parent),
objectiveValue(sourceNode.objectiveValue),
childrenAlive(sourceNode.childrenAlive),
nodeNumber(sourceNode.nodeNumber),
nodeDepth(sourceNode.nodeDepth) {}
BBSMPSNode::BBSMPSNode(BBSMPSNode* parent_ptr, std::vector<BBSMPSBranchingInfo>& bInfos){

  
  if (parent_ptr!=NULL){
    parent=parent_ptr;
    parent->incrementAliveChildren();
    nodeDepth=parent->getNodeDepth()+1;
  }
  else {
    parent=NULL;
    nodeDepth=0;

  }

  objectiveValue=-INFINITY;
  branchingInfos=std::vector<BBSMPSBranchingInfo>(bInfos);
  childrenAlive=0;
  nodeNumber=(++nodeCounter);
}

BBSMPSNode::~BBSMPSNode(){
  assert(childrenAlive==0);
  if(parent!=NULL){
    parent->decrementAliveChildren();
  }
  
}

double BBSMPSNode::getParentObjective() const{
	if (parent!=NULL){
    return parent->getObjective();
  }
  return -INFINITY;
}

double BBSMPSNode::getObjective() const{
  return objectiveValue;
}

void BBSMPSNode::setParentPtr(BBSMPSNode* parent_ptr){
  parent=&(*parent_ptr);
}

BBSMPSNode* BBSMPSNode::getParentPtr(){
  return parent;
}



void BBSMPSNode::setObjective(double lb){
  objectiveValue=lb;
}


void BBSMPSNode::addBranchingInformation(BBSMPSBranchingInfo& bi){
  branchingInfos.push_back(bi);
}


void BBSMPSNode::auxCopyAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector){
  if (parent!=NULL)parent->auxCopyAllBranchingInformation(biVector);
  if (branchingInfos.size()>0)  biVector.insert(biVector.end(),branchingInfos.begin(), branchingInfos.end());
  
}

void BBSMPSNode::incrementAliveChildren(){
  childrenAlive++;
} 

void BBSMPSNode::decrementAliveChildren(){
  assert(childrenAlive>0);//we should always have children before calling this.

  childrenAlive--;
  if (childrenAlive==0){

    delete this;
  }


}


void BBSMPSNode::eliminate(){

  if (childrenAlive==0){

    delete this;
  }

}

void BBSMPSNode::setIncrementalWarmStartState(std::vector< std::pair < BAIndex, variableState > > &changes){
  partialStartState=changes;
}


void BBSMPSNode::reconstructWarmStartState(BAFlagVector<variableState> &state){

  if (parent!=NULL) parent->reconstructWarmStartState(state);
  for (int i=0; i< partialStartState.size(); i++){
    //std::cout<<"Updating "<<partialStartState[i].first.scen<<" "<<partialStartState[i].first.idx<<" with "<<partialStartState[i].second<<endl;
    state.getVec(partialStartState[i].first.scen)[partialStartState[i].first.idx]=partialStartState[i].second;
  }
}




int BBSMPSNode::getBranchingInfoSize(){
  return branchingInfos.size();
}


void BBSMPSNode::getAllBranchingInformation(std::vector<BBSMPSBranchingInfo> &biVector){
  BBSMPSNode* n_ptr = parent;
  int branchingInfoSize=getBranchingInfoSize();
  while(n_ptr!=NULL){
    branchingInfoSize+=n_ptr->getBranchingInfoSize();
    n_ptr=n_ptr->parent;
  }
  biVector.reserve(branchingInfoSize);
  auxCopyAllBranchingInformation(biVector);
}

void BBSMPSNode::getAllBranchingInformation(denseBAVector &lb,denseBAVector &ub){

  std::vector<BBSMPSBranchingInfo> allBranchings;

  getAllBranchingInformation(allBranchings);

  for (int i=0; i<allBranchings.size(); i++){
    if (allBranchings[i].getStageNumber()==1){
     int varN=allBranchings[i].getVarNumber();
     double bd=allBranchings[i].getBound();
     if (allBranchings[i].getDirection()=='L'){
      lb.getFirstStageVec()[varN] = std::max(bd,lb.getFirstStageVec()[varN]);
    }
    else if (allBranchings[i].getDirection()=='U'){
      ub.getFirstStageVec()[varN] =  std::min(bd,ub.getFirstStageVec()[varN]);
    }
    else if(allBranchings[i].getDirection()=='E'){
      lb.getFirstStageVec()[varN] = std::max(bd,lb.getFirstStageVec()[varN]);
      ub.getFirstStageVec()[varN] =  std::min(bd,ub.getFirstStageVec()[varN]);
    }
  }
  else{
    int scenario=allBranchings[i].getScenarioNumber();
    int varN=allBranchings[i].getVarNumber();
    double bd=allBranchings[i].getBound();
    if (allBranchings[i].getDirection()=='L'){

      lb.getSecondStageVec(scenario)[varN] = std::max(bd,lb.getSecondStageVec(scenario)[varN]);
    }
    else if (allBranchings[i].getDirection()=='U'){
      ub.getSecondStageVec(scenario)[varN] =  std::min(bd,ub.getSecondStageVec(scenario)[varN]);
    }
    else if(allBranchings[i].getDirection()=='E'){
      lb.getSecondStageVec(scenario)[varN] = std::max(bd,lb.getSecondStageVec(scenario)[varN]);
      ub.getSecondStageVec(scenario)[varN] =  std::min(bd,ub.getSecondStageVec(scenario)[varN]);
    }
  }
}


}

int BBSMPSNode::getNodeNumber()const {
  return nodeNumber;
}

void  BBSMPSNode::setNodeDepth(int depth){
  nodeDepth=depth;
}

int BBSMPSNode::getNodeDepth(){
  return nodeDepth;
}





