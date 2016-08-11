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
#include "BBSMPSCuttingPlaneGeneratorGMI.hpp"

using namespace std;

void fillPlaneExpression(denseBAVector &expr, sparseBAVector &row, double beta){

	const BADimensionsSlacks &originalDims= BBSMPSSolver::instance()->getOriginalBADimensionsSlacks();

	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    double bHat= floorFracPart(beta);
	const CoinIndexedVector &v1 = row.getVec(-1).v;
	const double *v1Elts = v1.denseVector();
	const int* v1Idx = v1.getIndices();
	int nnz1 = v1.getNumElements();
	int numFirstStageVars = dimsSlacks.inner.numVars(-1);
	int originalFirstStageVars = originalDims.inner.numVars(-1);
	for (int j = 0; j < nnz1; j++) {
		int r = v1Idx[j];
		double val = v1Elts[r];
		if (val!=0){
			double coef=0;
			if (r<originalFirstStageVars && input.isFirstStageColInteger(r)){
				double aHat=floorFracPart(val);
				if (aHat<bHat)coef= aHat/bHat;
				else coef= ((1-aHat))/(1-bHat);
			}
			else {
				if (val>0) coef=val/bHat;
				else coef= (-val/(1-bHat));
			}
			if (coef!=0)expr.getFirstStageVec()[r]=coef;
		}

	}
	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			const CoinIndexedVector &v2 = row.getSecondStageVec(scen).v;
			const double *v2Elts = v2.denseVector();
			const int* v2Idx = v2.getIndices();
			int nnz2 = v2.getNumElements();
			int numSecStageVars = dimsSlacks.inner.numVars(scen);
			int originalSecondStageVars = originalDims.inner.numVars(scen);

			for (int j = 0; j < nnz2; j++) {
				int r = v2Idx[j];
				double val = v2Elts[r];
				if (val!=0){
					double coef=0;
					if (r<originalSecondStageVars &&input.isSecondStageColInteger(scen,r)){
						double aHat=floorFracPart(val);
						if (aHat<bHat)coef= aHat/bHat;
						else coef= ((1-aHat))/(1-bHat);
					}
					else {
						if (val>0) coef=val/bHat;
						else coef= (-val/(1-bHat));
					}
					if (coef!=0)expr.getSecondStageVec(scen)[r]=coef;
				}

			}
		}
	}



}

bool BBSMPSCuttingPlaneGeneratorGMI::generateCuttingPlane(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	//Generate cut here
	//BBSMPS_ALG_LOG_SEV(info) <<"We are generating GMI cuts"<< status;
	PIPSSInterface &solver= BBSMPSSolver::instance()->getPIPSInterface();
    const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
    BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
    SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();

    sparseBAVector beta;
	beta.allocate(dimsSlacks,ctx,BasicVector);
	beta.clear();
	solver.generateBetas(beta);

	int count=0;
	const CoinIndexedVector &v1 = beta.getFirstStageVec().v;
	const double *v1Elts = v1.denseVector();
	const int* v1Idx = v1.getIndices();
	int nnz1 = v1.getNumElements();

	for (int j = 0; j < nnz1; j++) {
		int row = v1Idx[j];
		double val = v1Elts[row];
		if (!isIntFeas(val, intTol)){
			sparseBAVector newRow;
			BAIndex inIndex;
			inIndex.scen=-1;
			inIndex.idx=row;
			solver.generateNonBasicRow(inIndex, newRow);
			count++;

			denseBAVector expr;
			expr.allocate(dimsSlacks,ctx,PrimalVector);
			expr.clear();

			fillPlaneExpression(expr, newRow, val);
			double bHat= floorFracPart(val);

			BBSMPSCuttingPlane *cuttingPlane= new BBSMPSCuttingPlane(bHat, COIN_DBL_MAX, expr);
			node->addCuttingPlane(cuttingPlane);
		}

	}


	for (int scen = 0; scen < input.nScenarios(); scen++) {
		if(ctx.assignedScenario(scen)) {
			const CoinIndexedVector &v2 = beta.getSecondStageVec(scen).v;
			const double *v2Elts = v2.denseVector();
			const int* v2Idx = v2.getIndices();
			int nnz2 = v2.getNumElements();
			for (int j = 0; j < nnz2; j++) {
				int row = v2Idx[j];
				double val = v2Elts[row];
				if (!isIntFeas(val, intTol)){
					sparseBAVector newRow;

					BAIndex inIndex;
					inIndex.scen=scen;
					inIndex.idx=row;
					solver.generateNonBasicRow(inIndex, newRow);
					count++;

					denseBAVector expr;
					expr.allocate(dimsSlacks,ctx,PrimalVector);
					expr.clear();

					fillPlaneExpression(expr, newRow, val);
					double bHat= floorFracPart(val);

					BBSMPSCuttingPlane *cuttingPlane= new BBSMPSCuttingPlane(1, COIN_DBL_MAX, expr);
					node->addCuttingPlane(cuttingPlane);

				}
			}
		}
	}



	return false;

}

bool BBSMPSCuttingPlaneGeneratorGMI::shouldItRun(BBSMPSNode* node, denseBAVector &LPRelaxationSolution){
	return (node->getNodeDepth()<2);

}
