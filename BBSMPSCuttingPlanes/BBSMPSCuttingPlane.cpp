#include "BBSMPSCuttingPlane.hpp"

using namespace std;

BBSMPSCuttingPlane::BBSMPSCuttingPlane(double _lb, double _ub, denseBAVector &_expr){
	lb=_lb;
	ub=_ub;
	expr=denseBAVector(_expr);
	//cout<<"Created a new cutting plane!!!!!! The size of he expression is "<<expr.getFirstStageVec().length()<<endl;
}
BBSMPSCuttingPlane::~BBSMPSCuttingPlane(){};

bool BBSMPSCuttingPlane::applyCuttingPlane(){




/*

	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	cout<<"ch1 "<<endl;
	//Create as many new vars as scenarios
	
 cout<<"ch2 "<<endl;
	//Create new vectors for first phase: copy the elements, set to ones the new variables
 	//int firstIndex=rootSolver.addFirstStageColumn(0,0,0);
//cout<<" index of new var "<<firstIndex<<endl;
 	int originalNumVars1= dimsSlacks.inner.numFirstStageVars();
	
 	vector<double> Acons(originalNumVars1,1);



 	/*Acons[0]=1;
 	Acons[1]=0;

 	Acons[3]=0;
 	Acons[4]=0;

 	Acons[firstIndex]=1;

 	//Tcons[1]=1;
 	Wcons[0]=1;

 	Tcons[0]=1;
 	
 	Tcons[2]=1;
 	
 	//Acons[5]=1;
 	
 	Tcons[firstIndex]=1;*/
	//At this point, the A vector is ready
 	//vector< vector <double> > elts1;
 	//elts1.push_back(Acons);
 	//vector<double> firstStageRowsLb(1,COIN_DBL_MIN);
 	//vector<double> firstStageRowsUb(1,COIN_DBL_MAX);
 	//rootSolver.addFirstStageRows(elts1, firstStageRowsLb, firstStageRowsUb, 1);
/*
 	int numScenarios=input.nScenarios();
 	for (int scen=0; scen<numScenarios; scen++){
	 	if(ctx.assignedScenario(scen)){
	 		int originalNumVars2= dimsSlacks.inner.numSecondStageVars(0);
	

 	vector<double> Wcons(originalNumVars2,1);
 	vector<double> Tcons(originalNumVars1,1);
rootSolver.addRow(Tcons,Wcons,scen,COIN_DBL_MIN,COIN_DBL_MAX);
 	

	 	}
	 }
 	rootSolver.addFirstStageRow(Acons,COIN_DBL_MIN,COIN_DBL_MAX);
 	//rootSolver.commitNewRows();


/*




	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int originalNumVars1= dimsSlacks.inner.numFirstStageVars();

	//Create as many new vars as scenarios
	int numScenarios=input.nScenarios();
 	int firstIndex=-1;
 	if(numScenarios>0){
 		firstIndex=rootSolver.addFirstStageColumn(0,COIN_DBL_MAX,0);
 		for(int i=1; i< numScenarios; i++)rootSolver.addFirstStageColumn(0,COIN_DBL_MAX,0);
 	}
	

	//Create new vectors for first phase: copy the elements, set to ones the new variables
 	cout<<originalNumVars1<<" and now "<<dimsSlacks.inner.numFirstStageVars()<<endl;
 	vector<double> Acons(dimsSlacks.inner.numFirstStageVars(),1);

	
	//At this point, the A vector is ready
 	//vector< vector <double> > elts1;
 	//elts1.push_back(Acons);
 	//vector<double> firstStageRowsLb(1,lb);
 	//vector<double> firstStageRowsUb(1,ub);
 	//rootSolver.addFirstStageRows(elts1, firstStageRowsLb, firstStageRowsUb, 1);
 	rootSolver.addFirstStageRow(Acons,lb, ub);//COIN_DBL_MIN,COIN_DBL_MAX); //lb,ub);
	//Create new vectors for each of the scenarios
	//For each scenario:
	for (int scen=0; scen<numScenarios; scen++){
	 	if(ctx.assignedScenario(scen)){

			//Generate T with ones on the variable that is local
			vector<double> Tcons(dimsSlacks.inner.numFirstStageVars(),0);
			Tcons[firstIndex+scen]=-1;
			//Generate W with the elements belonging to the second expression
			int numVars2= dimsSlacks.inner.numSecondStageVars(scen);
			
			vector<double> Wcons(numVars2,1);
			

			//Set bounds to zero and add new constraint.
			//vector< vector <double> > Telts2;
			//Telts2.push_back(Tcons);

			//vector< vector <double> > Welts2;
			//Welts2.push_back(Wcons);
		
			//vector<double> SecondStageRowsLb(1,0);
			//vector<double> SecondStageRowsUb(1,0);	
			rootSolver.addRow(Tcons,Wcons,scen,0,0);
			//rootSolver.addSecondStageRows(Telts2, Welts2, scen, SecondStageRowsLb, SecondStageRowsUb, 1);

	 		}


	 	}
*/




	SMPSInput &input =BBSMPSSolver::instance()->getSMPSInput();
	PIPSSInterface &rootSolver= BBSMPSSolver::instance()->getPIPSInterface();
	const BADimensionsSlacks &dimsSlacks= BBSMPSSolver::instance()->getBADimensionsSlacks();
	BAContext &ctx=BBSMPSSolver::instance()->getBAContext();
	int originalNumVars1= dimsSlacks.inner.numFirstStageVars();
	int originalTotalSize= dimsSlacks.numFirstStageVars();
	//Create as many new vars as scenarios
	int numScenarios=input.nScenarios();
 	int firstIndex=-1;
 	if(numScenarios>0){
 		firstIndex=rootSolver.addFirstStageColumn(0,1,1);
 		for(int i=1; i< numScenarios; i++)rootSolver.addFirstStageColumn(0,1,1);
 	}
	

	//Create new vectors for first phase: copy the elements, set to ones the new variables
 	cout<<originalNumVars1<<" and now "<<dimsSlacks.inner.numFirstStageVars()<<endl;
 	vector<double> Acons(dimsSlacks.inner.numFirstStageVars(),1);

 	const denseVector &cutFirstStageVec = expr.getFirstStageVec();
 	 cout<<"ch22 "<<cutFirstStageVec.length()<<" "<<originalNumVars1<<" "<<originalTotalSize<< endl;
 	for(int i=0; i< originalNumVars1; i++)Acons[i]=0;//cutFirstStageVec[i];

	cout<<"The expression looks like that"<<endl;
cout<<lb<<" "<<ub<< "Acol ";
	for (int i=0; i< Acons.size(); i++)cout<<"["<<i<<"]"<<Acons[i]<<" ";
		cout<<endl;


cout<<"ch3 "<<endl;
	//At this point, the A vector is ready
 	//vector< vector <double> > elts1;
 	//elts1.push_back(Acons);
 	//vector<double> firstStageRowsLb(1,lb);
 	//vector<double> firstStageRowsUb(1,ub);
 	//rootSolver.addFirstStageRows(elts1, firstStageRowsLb, firstStageRowsUb, 1);
 	rootSolver.addFirstStageRow(Acons,0,0);//lb,ub);// //
cout<<"ch4 "<<endl;
	//Create new vectors for each of the scenarios
	//For each scenario:
	for (int scen=0; scen<numScenarios; scen++){
	cout<<"ch5 "<<endl;
	 	if(ctx.assignedScenario(scen)){

			//Generate T with ones on the variable that is local
			vector<double> Tcons(dimsSlacks.inner.numFirstStageVars(),0);
		//	cout<<" setting var numver "<<firstIndex+scen<<" in T "<<endl;
			//Tcons[firstIndex+scen]=-1;

	//	/	cout<<"Tcol "<<scen<<" ";
	//for (int i=0; i< Tcons.size(); i++)cout<<"["<<i<<"]"<<Tcons[i]<<" ";
	//	cout<<endl;


			//Generate W with the elements belonging to the second expression
			int numVars2= dimsSlacks.inner.numSecondStageVars(scen);
		//	cout<<" Stage has "<<numVars2<<" vars"<<endl;
			
			vector<double> Wcons(numVars2,0);
			const denseVector &cutSecondStageVec = expr.getSecondStageVec(scen);
			for(int i=0; i< numVars2; i++)Wcons[i]=cutSecondStageVec[i];


	//			cout<<"Wcons "<<scen<<" ";
	//for (int i=0; i< Wcons.size(); i++)cout<<"["<<i<<"]"<<Wcons[i]<<" ";
	//	cout<<endl;

			//Set bounds to zero and add new constraint.
			//vector< vector <double> > Telts2;
			//Telts2.push_back(Tcons);

			//vector< vector <double> > Welts2;
			//Welts2.push_back(Wcons);
		
			//vector<double> SecondStageRowsLb(1,0);
			//vector<double> SecondStageRowsUb(1,0);	
			//rootSolver.addRow(Tcons,Wcons,scen,COIN_DBL_MIN,COIN_DBL_MAX);
			//rootSolver.addSecondStageRows(Telts2, Welts2, scen, SecondStageRowsLb, SecondStageRowsUb, 1);

	 		}

 	}
cout<<"ch6 "<<endl;
		
}