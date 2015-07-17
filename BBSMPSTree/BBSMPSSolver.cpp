#include "BBSMPSSolver.hpp"
using namespace std;

BBSMPSSolver *BBSMPSSolver::solverInstance = 0;


BAContext& BBSMPSSolver::getBAContext(){
  return ctx;
}

int& BBSMPSSolver::getMype(){
  return mype;
}

SMPSInput& BBSMPSSolver::getSMPSInput(){
  return input;
}

PIPSSInterface& BBSMPSSolver::getPIPSInterface(){
  return rootSolver;
}

BADimensions& BBSMPSSolver::getBADimensions(){
  return dims; 
}

BADimensionsSlacks& BBSMPSSolver::getBADimensionsSlacks(){
  return dimsSlacks;
}

denseBAVector& BBSMPSSolver::getOriginalLB(){
  return lb;
}

denseBAVector& BBSMPSSolver::getOriginalUB(){
  return ub;
}

BBSMPSSolver *BBSMPSSolver::instance(){
  return solverInstance;
}

//TODO: The relaxation could not be initialized. we need proper guards.
denseBAVector& BBSMPSSolver::getLPRelaxation(){
  return LPRelaxation;
}
void BBSMPSSolver::setLPRelaxation(denseBAVector &_LPRelaxation){
  LPRelaxation=_LPRelaxation;
}


BBSMPSSolver *BBSMPSSolver::initialize(const SMPSInput &_input){
 if (solverInstance) delete solverInstance;
 solverInstance = new BBSMPSSolver(_input);
 return solverInstance;
}

BBSMPSSolver::BBSMPSSolver(const SMPSInput &_input):
ctx(MPI_COMM_WORLD),
mype(ctx.mype()),
input(_input),
problemData(input,ctx),
rootSolver(problemData, PIPSSInterface::useDual),
pre(problemData,input),
dims(problemData.dims.inner),
dimsSlacks(dims),
startTimeStamp(MPI_Wtime()){
  lb = rootSolver.getLB();
  ub = rootSolver.getUB();

}

double BBSMPSSolver::getWallTime(){
  return MPI_Wtime()-startTimeStamp;
}

bool BBSMPSSolver::isInitialized(){
  return (solverInstance!=NULL);
}

void BBSMPSSolver::deInitialize(){
  if(isInitialized()){
    delete solverInstance;
  }
}

const BAFlagVector<variableState>& BBSMPSSolver::getOriginalWarmStart(){
  return originalWarmStart;
}
void BBSMPSSolver::setOriginalWarmStart(BAFlagVector<variableState>&warmStart){
  originalWarmStart=warmStart;
}

void BBSMPSSolver::printPresolveStatistics(){


  BBSMPS_ALG_LOG_SEV(warning)<<"~~~~~~~~~~~~~~~PRESOLVE STATISTICS~~~~~~~~~~~~~~~~";
   
    BBSMPS_ALG_LOG_SEV(warning)<<"Number Of Presolves:"<<pre.numPresolves<<":Number Of Bound Chgs:"<<pre.numBdsChg;
  BBSMPS_ALG_LOG_SEV(warning)<<"Number Of RHS Chgs:"<<pre.numRhsChg<<":Number Of Coefficient Chgs:"<<pre.numCoeffChg;
   

BBSMPS_ALG_LOG_SEV(warning)<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
}


void BBSMPSSolver::addSolutionToPool(BBSMPSSolution &sol){
    solutionPool.insert(sol);
};


void BBSMPSSolver::printSolutionStatistics(double objLB){


    double bestTime=COIN_DBL_MAX;

    for (std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin(); it!=solutionPool.end(); ++it){
      BBSMPSSolution s = *it;
      if (s.getTimeOfDiscovery()<bestTime)bestTime=s.getTimeOfDiscovery();
    }

      BBSMPS_ALG_LOG_SEV(warning)<<"---------------SOLUTION STATISTICS----------------";
      
      BBSMPS_ALG_LOG_SEV(warning)<<"Solution Pool Size:"<<solutionPool.size()<<":Time To First Solution:"<<bestTime;
  int itCounter=0;
  for (std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin(); it!=solutionPool.end(); ++it){
      BBSMPSSolution s = *it;
      double solGap = fabs(s.getObjValue()-objLB)*100/(fabs(objLB)+10e-10);
      BBSMPS_ALG_LOG_SEV(warning)<<"Solution:"<<itCounter<<":Solution Value:"<<s.getObjValue()<<":Time Of Discovery:"<<s.getTimeOfDiscovery()<<":Solution Gap:"<<solGap;
    itCounter++;
    }
  
  BBSMPS_ALG_LOG_SEV(warning)<<"--------------------------------------------------";
}

const BBSMPSSolution &BBSMPSSolver::getSoln(int index){
  assert(index<solutionPool.size());

  std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin();
  int ctr=0;
  while (ctr<index && ctr< solutionPool.size()){
    it++;
    ctr++;
  }
  return (*it);

}
const BBSMPSSolution &BBSMPSSolver::getSolnBySolNumber(int number){
  
  std::multiset<BBSMPSSolution,solutionComparison>::iterator it=solutionPool.begin();
  int ctr=0;
  while ( ctr< solutionPool.size()){
    if ((*it).getSolNumber()==number)return (*it);
    it++;
    ctr++;
  }
  assert(false);
  return (*solutionPool.begin());

}



int BBSMPSSolver::getSolPoolSize(){
  return solutionPool.size();
}



