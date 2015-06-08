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
rootSolver(input, ctx, PIPSSInterface::useDual),
dims(input, ctx),
dimsSlacks(dims){
  lb = rootSolver.getLB();
  ub = rootSolver.getUB();
}

bool BBSMPSSolver::isInitialized(){
  return (solverInstance!=NULL);
}


