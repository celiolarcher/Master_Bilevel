#include "LagrangeMultpSimplex.h"

#include <cmath>
#include <iostream>
using namespace std;

#define TOL_ACT_CONST 1e-2

extern double TOL_EQ_CONST;
extern double TOL_NEQ_CONST;

extern bool InfeasibleAvaliation;

    int LagrangeMultpSimplex::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW()+1;

	constraintsNEQNumber=function->getNEQConstraintNumberUP()+function->getNEQConstraintNumberLW();
	constraintsEQNumber=constraintsNumber-constraintsNEQNumber;
	editBegin=0;
	editSize=function->getDimensionUP();

	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}

	return 1;
    }

    int LagrangeMultpSimplex::decodifySolution(Solution *sol){
      
          sol->feasible=1;
	  sol->completeSolution=1;
          
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

          if(sol->feasible)// || InfeasibleAvaliation)
	 sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
          
          return 1;
    }

#include <iostream>
using namespace std;
     int LagrangeMultpSimplex::simplexLagrangeMultipliers(double x[], double y[], double h[],double *opt){
		int mark[function->getNEQConstraintNumberLW()];
		int countMultp=0;
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		        if(h[i]<-TOL_ACT_CONST || h[i]>TOL_ACT_CONST) mark[i]=0;
			else{ 	
				mark[i]=1;
				countMultp++;
			}
		}


		int sizeCol=(function->getNEQConstraintNumberLW()+2*function->getDimensionLW()+1);

		double tableau[(function->getDimensionLW()+1)*sizeCol];
		function->getSimplexTableauKKT(x,y,tableau);		

		int endCol=countMultp+2*function->getDimensionLW()+1;
/*
		for(int i=0;i<function->getDimensionLW()+1;i++){
			for(int j=0;j<sizeCol;j++) cout<<tableau[i*sizeCol+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
*/
		int swapPoint=(function->getNEQConstraintNumberLW()+1)*function->getDimensionLW()-1;
		for(int i=function->getDimensionLW()-1;i>=0+1;i--){
			for(int j=function->getNEQConstraintNumberLW()+1-1;j>=0;j--) tableau[i*sizeCol+j]=tableau[swapPoint--];

		}

		
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
			for(int j=0;j<endCol;j++) cout<<tableau[i*sizeCol+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
*/		

		int loop=1;

//if(countMultp==0)cout<<"+";	
		while(loop++){
			double valueRef=-1;
			int pivot=-1;
			for(int j=0;j<endCol-1;j++){
				if(tableau[function->getDimensionLW()*sizeCol+j]>valueRef){
					valueRef=tableau[function->getDimensionLW()*sizeCol+j];
					pivot=j;
				}
			}

			if(valueRef<=1e-10) break;  //tolerancia

			int line=-1;
			valueRef=-1;
			for(int i=0;i<function->getDimensionLW();i++){
				if(tableau[i*sizeCol+pivot]>0 && (valueRef<0 || valueRef>(tableau[i*sizeCol+endCol-1]/tableau[i*sizeCol+pivot]))){
					line=i;
					valueRef=tableau[i*sizeCol+endCol-1]/tableau[i*sizeCol+pivot];
				}
			}

			if(valueRef<0) {
				tableau[function->getDimensionLW()*sizeCol+endCol-1]=1e5;cout<<"erro";
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
				tableau[function->getDimensionLW()*sizeCol+endCol-1]=1e5;cout<<"erro";
				break;
			}
		}

		(*opt)=tableau[function->getDimensionLW()*sizeCol+endCol-1];
/*
		//if((*opt)==0){
			for(int i=0;i<function->getDimensionLW()+1;i++){
					for(int j=0;j<endCol;j++) cout<<tableau[i*sizeCol+j]<<"\t";
					cout<<"\n";
				}
			*/	
		
/*
		double wVec[function->getDimensionLW()];
		
		for(int j=endCol-1-2*function->getDimensionLW(), wCount=0;j<endCol-1-function->getDimensionLW();j++,wCount++){
		    int flag=-2;
		  
		    wVec[wCount]=0;
		    for(int i=0;i<function->getDimensionLW()+1;i++){
		        if(tableau[i*sizeCol+j]!=0 && flag!=-2) flag=-1;
		        if(tableau[i*sizeCol+j]!=0 && flag==-2) flag=i;
		    }
		    
		    if(flag>=0) wVec[wCount]=tableau[flag*sizeCol+endCol-1];
		    else{
		        flag=-2;
		        for(int i=0;i<function->getDimensionLW()+1;i++){
		          if(tableau[i*sizeCol+j+function->getDimensionLW()]!=0 && flag!=-2) flag=-1;
		          if(tableau[i*sizeCol+j+function->getDimensionLW()]!=0 && flag==-2) flag=i;
		        }
		        if(flag>=0) wVec[wCount]=tableau[flag*sizeCol+endCol-1];	
		    }
		}
		     */
	
		 /*
		double aux=0;
		for(int i=0;i<function->getDimensionLW();i++) aux+=wVec[i];*/
		/*
		 
		if(aux+1e-10<(*opt) || aux-1e-10>(*opt)){
		  
		  			for(int i=0;i<function->getDimensionLW()+1;i++){
					for(int j=0;j<endCol;j++) cout<<tableau[i*sizeCol+j]<<"\t";
					cout<<"\n";
				}
		    cout<<aux<<"\t"<<(*opt)<<"\n";
		 
		}*/
		/*
		aux=0;
		for(int i=0;i<function->getDimensionLW();i++) aux+=wVec[i]*wVec[i];
		
		(*opt)=sqrt(aux);*/
	//	}else{cout<<(*opt);}
	/*
		aux=wVec[0];
		for(int i=1;i<function->getDimensionLW();i++) {
		  if(aux<wVec[i])
		      aux=wVec[i];
		}
		(*opt)=aux;
	        */
		return 1;
     }


