#include "DELowerLevel.h"
#include "DEBPenalty.h"

#include <cmath>
#include <iostream>
using namespace std;

#define TOL_LEQ_CONST 1e-5
#define PENALTY_DEF 1e5

extern double TOL_EQ_CONST;
extern double TOL_NEQ_CONST;

extern bool InfeasibleAvaliation;

    int DELowerLevel::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();

	constraintsNEQNumber=function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW(); //Lemke é executado no início do teste de restrições e considerado de não igualdade
	constraintsEQNumber=constraintsNumber-constraintsNEQNumber;
	editBegin=0;
	editSize=function->getDimensionUP();
	
	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}
	
	UP=this;
	
	LW=new DELowerLevel();
	LW->function=function;
	
	LW->solutionSize=function->getDimensionUP()+function->getDimensionLW();
	LW->constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();

	LW->constraintsNEQNumber=function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW(); //Lemke é executado no início do teste de restrições e considerado de não igualdade
	LW->constraintsEQNumber=constraintsNumber-constraintsNEQNumber;
	LW->editBegin=function->getDimensionUP();
	LW->editSize=function->getDimensionLW();
	
	LW->boundAttributes=boundAttributes;
	
	DELW=NULL;


	return 1;
    }

    int DELowerLevel::decodifySolution(Solution *sol){
      
          sol->feasible=1;
          sol->completeSolution=1;
          
          sol->upLevelFunction=0;
	  
          if(sol->countConstraint==0) return 1;
          
          int offset=0;
        
          if(this==this->UP){
	PenaltySolution *penalty=new DEBPenalty();
	
	if(!DELW)
	  DELW=new DifferentialEvolution(this->LW,penalty,30);  
	else
	 DELW->resetPopulation();
	
	for(int i=0;i<30;i++){
	    for(int j=this->UP->editBegin;j<this->UP->editBegin+this->UP->editSize;j++){
	        DELW->Population[i]->vectorCharacters[j]=sol->vectorCharacters[j];
	        DELW->nextPopulation[i]->vectorCharacters[j]=sol->vectorCharacters[j];
	    }
	}
	
	for(int i=0;i<100;i++){
	    for(int p=0;p<30;p++)
	        DELW->mutatePopulation_TargetToRand_1_Wall(0.7,0.7,p);
	    
	    for(int p=0;p<30;p++)
	        DELW->recombinePopulation(0.9,p);
	    
	    for(int p=0;p<30;p++)
	        DELW->selectPopulation(p);
	}
	for(int j=this->LW->editBegin;j<this->LW->editBegin+this->LW->editSize;j++)sol->vectorCharacters[j]=0;

	
	if(DELW->best && DELW->best->feasible)
	  for(int j=this->LW->editBegin;j<this->LW->editBegin+this->LW->editSize;j++) sol->vectorCharacters[j]=DELW->best->vectorCharacters[j];
	else{
	        sol->completeSolution=0;
	        sol->feasible=0;
	        for(int i=0;i<function->getDimensionUP();i++)
	          cout<<sol->vectorCharacters[i]<<"\t";
	          
	        cout<<"\n";
	        //for(int i=0;i<function->getDimensionLW();i++)
	         // sol->vectorCharacters[i+function->getDimensionUP()]=(function->bounds[2*(i+function->getDimensionUP())]+function->bounds[2*(i+function->getDimensionUP())+1])/2.0;//avg[i+function->getDimensionUP()];
	}
	  
          //DELW->clearPopulation();
          //delete DELW;
          //offset++;

	//cout<<sol->constraintValues[offset-1]<<"\t";
          //if(sol->constraintValues[offset-1]>TOL_EQ_CONST) sol->feasible=0;
          
	  if(sol->completeSolution){
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

		//  if(sol->feasible)// || InfeasibleAvaliation)
		 	sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
	  }else{
		for(int i=0;i<constraintsNumber;i++)sol->constraintValues[i]=0;
		
		if(constraintsNumber>0){
		    for(int i=0;i<function->getDimensionUP();i++){
		        if(sol->vectorCharacters[i]<function->bounds[2*i])
			sol->constraintValues[0]+=function->bounds[2*i]-sol->vectorCharacters[i];
		        
		        if(sol->vectorCharacters[i]>function->bounds[2*i+1])
			sol->constraintValues[0]+=sol->vectorCharacters[i]-function->bounds[2*i+1];
		    }
		}
	  }

	  return 1;
          }else{
	      if(function->getNEQConstraintNumberLW()>0){
	        function->constraintsValueNEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
	        for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){
		      if(sol->constraintValues[i]>TOL_NEQ_CONST){
			    sol->feasible=0;
		      }
	        }
	      }
	      
	      offset+=function->getNEQConstraintNumberLW();
	      
	      if(function->getEQConstraintNumberLW()>0){
		      function->constraintsValueEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		      for(int i=offset;i<offset+function->getEQConstraintNumberLW();i++){
			    if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
			      sol->feasible=0;
			    }
		      }
	      }
	      
	      offset+=function->getEQConstraintNumberLW();

	//      if(sol->feasible)// || InfeasibleAvaliation)
		    sol->upLevelFunction=function->getLWLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
	      
	      return 1;
	
          }
    }