#include "LagrangeMultp.h"

#include <cmath>

    int LagrangeMultp::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW()+function->getKKTConstraintNumber()+function->getNEQConstraintNumberLW();

	constraintsNEQNumber=function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW();
	constraintsEQNumber=constraintsNumber-constraintsNEQNumber;


	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}
	
	for(int i=function->getDimensionUP()+function->getDimensionLW();i<solutionSize;i++){
	    boundAttributes[2*i]=0;
	    boundAttributes[2*i+1]=20;
	}	


	return 1;
    }

    int LagrangeMultp::decodifySolution(Solution *sol){

          sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
      
          sol->feasible=1;
          
          if(sol->countConstraint==0) return 1;
          
          int offset=0;

          double lambda[function->getNEQConstraintNumberLW()];
          
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
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
			sol->feasible=0;
		        }
		     //   if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
			lambda[i-offset]=sol->vectorCharacters[function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+(i-offset)];
		//        }else{
		//	lambda[i-offset]=0;
		  //      }
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
		  function->constraintsValueKKT(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW(),lambda,sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getKKTConstraintNumber();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getKKTConstraintNumber();

          if(function->getNEQConstraintNumberLW()>0){  //usar função q considera restrições ja calculadas
		  function->constraintsSlackness(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }

        
	return 1;
    }
