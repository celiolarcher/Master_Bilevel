#include "LagrangeMultpAPMSimplex.h"
#define TOL_EQ_CONST 10e-8
#define TOL_NEQ_CONST 10e-8


#include <cmath>
#include <iostream>
using namespace std;

    int LagrangeMultpAPMSimplex::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW()+1;

	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}
	/*
	for(int i=function->getDimensionUP()+function->getDimensionLW();i<solutionSize;i++){
	    boundAttributes[2*i]=0;
	    boundAttributes[2*i+1]=20;
	}	
*/

	return 1;
    }

    int LagrangeMultpAPMSimplex::decodifySolution(Solution *sol){

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
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }

          simplexLagrangeMultipliers(sol->vectorCharacters, sol->vectorCharacters+function->getDimensionUP(), sol->constraintValues+offset, sol->constraintValues+offset+function->getNEQConstraintNumberLW());
          
          offset+=function->getNEQConstraintNumberLW()+1;

	//cout<<sol->constraintValues[offset-1]<<"\t";
	  if(sol->constraintValues[offset-1]>TOL_EQ_CONST) sol->feasible=0;
          
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

	  return 1;
    }

    int LagrangeMultpAPMSimplex::updatePenalty(Solution **population, Solution **newPopulation, int sizePop, int sizeNewPop){
          levelUPMean=0;
          
          double violationMean[population[0]->countConstraint];
          if(!kWeight)
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
          
          if(newPopulation){
	  for(int i=0;i<sizeNewPop;i++){
	      levelUPMean+=newPopulation[i]->upLevelFunction;
	      for(int j=0;j<newPopulation[0]->countConstraint;j++){
		if(j < function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW()){
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
          }
          
          
          return 1;
     }

     int LagrangeMultpAPMSimplex::compareSolutions(Solution *sol1, Solution *sol2){   //1 se sol1 melhor que sol2, 0 caso contrário
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
	  //  if(sol1->upLevelFunction>levelUPMean)sol1->score+=sol1->upLevelFunction;
	  //  else sol1->score+=levelUPMean;
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
	 //   if(sol2->upLevelFunction>levelUPMean)sol2->score+=sol2->upLevelFunction;
	   // else sol2->score+=levelUPMean;
          }
	
          if(sol1->score < sol2->score) return 1;
          
          return 0;
     }


#include <iostream>
using namespace std;
     int LagrangeMultpAPMSimplex::simplexLagrangeMultipliers(double x[], double y[], double h[],double *opt){
		int mark[function->getNEQConstraintNumberLW()];
		int countMultp=0;
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		        if(h[i]<-TOL_EQ_CONST || h[i]>TOL_EQ_CONST) mark[i]=0;
			else{ 	
				mark[i]=1;
				countMultp++;
			}
		}


		int sizeCol=(function->getNEQConstraintNumberLW()+2*function->getDimensionLW()+1);

		double tableau[(function->getDimensionLW()+1)*sizeCol];
		function->getSimplexTableauKKT(x,y,tableau);		

		int endCol=countMultp+2*function->getDimensionLW()+1;


		for(int posMark=0,j=0;j<countMultp;j++,posMark++){
			
			for(;!mark[posMark];posMark++);
			for(int i=0;i<function->getDimensionLW();i++){
				tableau[i*sizeCol+j]=tableau[i*sizeCol+posMark];
			}
		}


		for(int i=0;i<function->getDimensionLW();i++){
			tableau[i*sizeCol+endCol-1]=tableau[i*sizeCol+function->getNEQConstraintNumberLW()];
		}

		for(int i=0;i<function->getDimensionLW();i++){
			for(int j=0;j<function->getDimensionLW();j++){
				tableau[i*sizeCol+j+countMultp]=(i==j);
			}
		}

		for(int i=0;i<function->getDimensionLW();i++){
			for(int j=0;j<function->getDimensionLW();j++){
				tableau[i*sizeCol+j+countMultp+function->getDimensionLW()]=-(i==j);
			}
		}

		for(int i=0;i<function->getDimensionLW();i++){
			if(tableau[i*sizeCol+endCol-1]<0){
				for(int j=0;j<endCol;j++) tableau[i*sizeCol+j]*=-1;
			}
		}
		
		for(int j=0;j<countMultp;j++)tableau[function->getDimensionLW()*sizeCol+j]=0;

		for(int j=countMultp;j<endCol-1;j++)tableau[function->getDimensionLW()*sizeCol+j]= -1;
		
		tableau[function->getDimensionLW()*sizeCol+endCol-1]=0;

		for(int i=0;i<function->getDimensionLW();i++){
			for(int j=0;j<endCol;j++)tableau[function->getDimensionLW()*sizeCol+j]+=tableau[i*sizeCol+j];
		}
		
/*
		for(int i=0;i<function->getDimensionLW()+1;i++){
			for(int j=0;j<endCol;j++) cout<<tableau[i][j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
*/
		int loop=1;
		while(loop++){
			double valueRef=-1;
			int pivot=-1;
			for(int j=0;j<endCol-1;j++){
				if(tableau[function->getDimensionLW()*sizeCol+j]>valueRef){
					valueRef=tableau[function->getDimensionLW()*sizeCol+j];
					pivot=j;
				}
			}

			if(valueRef<=0) break;

			int line=-1;
			valueRef=-1;
			for(int i=0;i<function->getDimensionLW();i++){
				if(tableau[i*sizeCol+pivot]>0 && (valueRef<0 || valueRef>(tableau[i*sizeCol+endCol-1]/tableau[i*sizeCol+pivot]))){
					line=i;
					valueRef=tableau[i*sizeCol+endCol-1]/tableau[i*sizeCol+pivot];
				}
			}

			if(valueRef<0) {
				tableau[function->getDimensionLW()*sizeCol+endCol-1]=10e5;
				break;
			}

			
			double multp=tableau[line*sizeCol+pivot];
			for(int j=0;j<endCol;j++){
				tableau[line*sizeCol+j]/=multp;
			}

			for(int i=0;i<function->getDimensionLW()+1;i++){
				if(i!=line){
					multp=tableau[i*sizeCol+pivot];
					for(int j=0;j<endCol;j++){
						tableau[i*sizeCol+j]-=multp*tableau[line*sizeCol+j];
					}
				}
			}

			if(loop>20) {
				tableau[function->getDimensionLW()*sizeCol+endCol-1]=10e5;cout<<"erro";
				break;
			}
		}

		(*opt)=tableau[function->getDimensionLW()*sizeCol+endCol-1];
/*
		if((*opt)==0){
			for(int i=0;i<function->getDimensionLW()+1;i++){
					for(int j=0;j<endCol;j++) cout<<tableau[i][j]<<"\t";
					cout<<"\n";
				}
		}else{cout<<(*opt);}
*/
		return 1;
     }


