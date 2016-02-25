#include "LagrangeMultpAPMSmooth.h"
#define TOL_EQ_CONST 10e-10
#define TOL_NEQ_CONST 10e-10
#define EPISILON 0




#include <cmath>
    double getLambda(Solution *sol, int pos, InputFunction *function){
	return sol->vectorCharacters[function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+pos];
    }

    int LagrangeMultpAPMSmooth::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW()+function->getKKTConstraintNumber();

	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}
	
	for(int i=function->getDimensionUP()+function->getDimensionLW();i<solutionSize;i++){
	    boundAttributes[2*i]=0;
	    boundAttributes[2*i+1]=10e5;
	}	


	return 1;
    }

    int LagrangeMultpAPMSmooth::decodifySolution(Solution *sol){

          sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
      
          sol->feasible=1;
          
          if(sol->countConstraint==0) return 1;
          
          int offset=0;

          if(function->getNEQConstraintNumberUP()>0){
		  function->constraintsValueNEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getNEQConstraintNumberUP();
          
          if(function->getNEQConstraintNumberLW()>0){
		  function->constraintsValueNEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){	
	  
			int multpLam=i-function->getNEQConstraintNumberUP();	
		
		  	sol->constraintValues[i]=(getLambda(sol,multpLam,function)-sol->constraintValues[i]) - sqrt((getLambda(sol,multpLam,function)+sol->constraintValues[i])*(getLambda(sol,multpLam,function)+sol->constraintValues[i])+4*EPISILON*EPISILON);

		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		               sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getNEQConstraintNumberLW();
          
          if(function->getEQConstraintNumberUP()>0){
		  function->constraintsValueEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
			  sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberUP();
          
          
          if(function->getEQConstraintNumberLW()>0){
		  function->constraintsValueEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberLW();
          
          if(function->getKKTConstraintNumber()>0){
		  function->constraintsValueKKT(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getKKTConstraintNumber();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
        
  	  return 1;
    }

#include <iostream>
using namespace std;
    int LagrangeMultpAPMSmooth::updatePenalty(Solution **population, Solution **newPopulation, int sizePop, int sizeNewPop){
          levelUPMean=0;
          
          double violationMean[population[0]->countConstraint];
          if(!kWeight)
		kWeight=new double[population[0]->countConstraint];
          
          for(int j=0;j<population[0]->countConstraint;violationMean[j]=0, j++);
	  
          for(int i=0;i<sizePop;i++){
		  levelUPMean+=population[i]->upLevelFunction;
		  for(int j=0;j<population[0]->countConstraint;j++){
			if(j < function->getNEQConstraintNumberUP()){
				if(population[i]->constraintValues[j]>0){
				    violationMean[j]+=population[i]->constraintValues[j];
				} 
			}else{
			 	violationMean[j]+=fabs(population[i]->constraintValues[j]);
			}
		  }
          }
          
          if(newPopulation){
	  for(int i=0;i<sizePop;i++){
		      levelUPMean+=newPopulation[i]->upLevelFunction;
		      for(int j=0;j<newPopulation[0]->countConstraint;j++){
			    if(j < function->getNEQConstraintNumberUP()){
				    if(newPopulation[i]->constraintValues[j]>0){
				        violationMean[j]+=newPopulation[i]->constraintValues[j];
				    } 
			    }else{
				    violationMean[j]+=fabs(newPopulation[i]->constraintValues[j]);
			    }
		      }
	  }
          }
          
          
          double violationSumSquare=0;
          for(int j=0;j<population[0]->countConstraint;j++){
		 violationMean[j]/=(2*sizePop);
		   
		 violationSumSquare+=violationMean[j]*violationMean[j];
          }
          
          levelUPMean/=(2*sizePop);
          
          for(int j=0;j<population[0]->countConstraint;j++){
	  	 kWeight[j]=(violationMean[j]*fabs(levelUPMean))/violationSumSquare;
	//	cout<<kWeight[j]<<"\t";
          }
          //cout<<endl;
     //     cout<<levelUPMean<<"\n";
          return 1;
     }

     int LagrangeMultpAPMSmooth::compareSolutions(Solution *sol1, Solution *sol2){   //1 se sol1 melhor que sol2, 0 caso contrário
         if(sol1->feasible) sol1->score=sol1->upLevelFunction;
          else{
	    sol1->score=0;
	    for(int j=0;j<sol1->countConstraint;j++){
	          if(j < function->getNEQConstraintNumberUP()){
			  if(sol1->constraintValues[j]>0)
			      sol1->score+=kWeight[j]*sol1->constraintValues[j];
	          }else{
		      sol1->score+=kWeight[j]*fabs(sol1->constraintValues[j]);
	          }
	    }

	    if(sol1->upLevelFunction>levelUPMean)sol1->score+=sol1->upLevelFunction;
	    else sol1->score+=levelUPMean;
          }
          
      
          if(sol2->feasible) sol2->score=sol2->upLevelFunction;
          else{
	    sol2->score=0;
	    for(int j=0;j<sol2->countConstraint;j++){
	          if(j < function->getNEQConstraintNumberUP()){
			  if(sol2->constraintValues[j]>0)
			      sol2->score+=kWeight[j]*sol2->constraintValues[j];
	          }else{
		      sol2->score+=kWeight[j]*fabs(sol2->constraintValues[j]);
	          }
	    }
	    if(sol2->upLevelFunction>levelUPMean)sol2->score+=sol2->upLevelFunction;
	    else sol2->score+=levelUPMean;
          }
	
          if(sol1->score < sol2->score) return 1;
          
          return 0;
     }



