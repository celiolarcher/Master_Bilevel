#include "OptLangrageMultp.h"
#define TOL_EQ_CONST 10e-5
#define TOL_NEQ_CONST 10e-5
#define APM

#include <iostream>
#include <cmath>
using namespace std;
    double levelUPMean;
    double *kWeight=NULL;

    double validateSolution(Solution *sol,InputFunction *function){  //Restrições NEQ seguidas das EQ
          sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
      
          sol->feasible=1;
          
          if(sol->countConstraint==0) return sol->feasible;
          
          int offset=0;

          if(function->getNEQConstraintNumberUP()>0){
		  function->constraintsValueNEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
		          sol->feasible=0;
		          #ifndef APM
			    return sol->feasible;
		          #endif
		        }
		  }
          }
          
          offset+=function->getNEQConstraintNumberUP();
          
          if(function->getNEQConstraintNumberLW()>0){
		  function->constraintsValueNEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
		          sol->feasible=0;
		          #ifndef APM
			    return sol->feasible;
		          #endif
		        }
		  }
          }
          
          offset+=function->getNEQConstraintNumberLW();
          
          if(function->getEQConstraintNumberUP()>0){
		  function->constraintsValueEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
			sol->feasible=0;
			#ifndef APM
			  return sol->feasible;
			#endif
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberUP();
          
          
          if(function->getEQConstraintNumberLW()>0){
		  function->constraintsValueEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		          #ifndef APM
			    return sol->feasible;
		          #endif
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberLW();
          
          if(function->getKKTConstraintNumber()>0){
		  function->constraintsValueKKT(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getKKTConstraintNumber();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		          #ifndef APM
			    return sol->feasible;
		          #endif
		        }
		  }
          }
          
          offset+=function->getKKTConstraintNumber();

          if(function->getNEQConstraintNumberLW()>0){  //usar função q considera restrições ja calculadas
		  function->constraintsSlackness(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		          #ifndef APM
			    return sol->feasible;
		          #endif
		        }
		  }
          }

          return sol->feasible;
     }


     int compareSolutions(Solution *sol1, Solution *sol2, InputFunction *function){  //1 se sol1 melhor que sol2, 0 caso contrário
          if((!sol1->feasible && sol2->feasible) || ((sol1->feasible && sol2->feasible || !sol1->feasible) && sol1->calcScore(function) >= sol2->calcScore(function)))
	  return 0;
          return 1;
     } 


     int updateAPM(Solution **population, int sizePop, InputFunction *function){
          levelUPMean=0;
          
          double violationMean[population[0]->countConstraint];
          kWeight=new double[population[0]->countConstraint];
          
          for(int j=0;j<population[0]->countConstraint;violationMean[j]=0, j++);
	  
          for(int i=0;i<sizePop;i++){
	  levelUPMean+=population[i]->upLevelFunction;
	  for(int j=0;j<population[0]->countConstraint;j++){
	        if(j < function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW()){
		if(population[i]->constraintValues[j]>0){
		    violationMean[j]+=population[i]->constraintValues[j];
		} 
	        }else{
		 violationMean[j]+=fabs(population[i]->constraintValues[j]);
	        }
	  }
          }
          
          
          double violationSumSquare=0;
          for(int j=0;j<population[0]->countConstraint;j++){
	   violationMean[j]/=sizePop;
	   
	   violationSumSquare+=violationMean[j]*violationMean[j];
          }
          
          levelUPMean/=sizePop;
          
          for(int j=0;j<population[0]->countConstraint;j++){
	   kWeight[j]=(violationMean[j]*levelUPMean)/violationSumSquare;
          }
          
          
          return 1;
     }

     int compareAPM(Solution *sol1, Solution *sol2, InputFunction *function){   //1 se sol1 melhor que sol2, 0 caso contrário
          if(sol1->feasible) sol1->score=sol1->upLevelFunction;
          else{
	    sol1->score=0;
	    for(int j=0;j<sol1->countConstraint;j++){
	          if(j < function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW()){
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
	          if(j < function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW()){
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
