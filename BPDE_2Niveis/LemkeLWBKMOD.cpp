#include "LemkeLW.h"

#include <cmath>
#include <iostream>
using namespace std;

#define TOL_LEQ_CONST 1e-5
#define PENALTY_DEF 1e5

extern double TOL_EQ_CONST;
extern double TOL_NEQ_CONST;

extern bool InfeasibleAvaliation;

    int LemkeLW::initInstance(InputFunction *function){
	this->function=function;

	solutionSize=function->getDimensionUP()+function->getDimensionLW();
	constraintsNumber=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+1;

	constraintsNEQNumber=function->getNEQConstraintNumberUP()+ 1; //Lemke é executado no início do teste de restrições e considerado de não igualdade
	constraintsEQNumber=constraintsNumber-constraintsNEQNumber;
	editBegin=0;
	editSize=function->getDimensionUP();
	
	boundAttributes=new double[2*solutionSize];
	for(int i=0;i<function->getDimensionUP();i++){
	    boundAttributes[2*i]=function->bounds[2*i];
	    boundAttributes[2*i+1]=function->bounds[2*i+1];
	}

	return 1;
    }

    int LemkeLW::decodifySolution(Solution *sol){
      
          sol->feasible=1;
          
          if(sol->countConstraint==0) return 1;
          
          int offset=0;
        
          getYLemke(sol->vectorCharacters, sol->vectorCharacters+function->getDimensionUP(), sol->constraintValues+offset);
          
          offset++;

	//cout<<sol->constraintValues[offset-1]<<"\t";
          if(sol->constraintValues[offset-1]>TOL_EQ_CONST) sol->feasible=0;
          
          if(function->getNEQConstraintNumberUP()>0){
		  function->constraintsValueNEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getNEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]>TOL_NEQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getNEQConstraintNumberUP();
      /*    
          if(function->getNEQConstraintNumberLW()>0){
	    function->constraintsValueNEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
	    for(int i=offset;i<offset+function->getNEQConstraintNumberLW();i++){
	          if(sol->constraintValues[i]>TOL_NEQ_CONST){
		sol->feasible=0;
	          }
	    }
          }
          
          offset+=function->getNEQConstraintNumberLW();
          */
          if(function->getEQConstraintNumberUP()>0){
		  function->constraintsValueEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberUP();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
			  sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberUP();
          /*
          if(function->getEQConstraintNumberLW()>0){
		  function->constraintsValueEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->constraintValues+offset);
		  for(int i=offset;i<offset+function->getEQConstraintNumberLW();i++){
		        if(sol->constraintValues[i]<-TOL_EQ_CONST || sol->constraintValues[i]>TOL_EQ_CONST){
		          sol->feasible=0;
		        }
		  }
          }
          
          offset+=function->getEQConstraintNumberLW();
*/
          if(sol->feasible)// || InfeasibleAvaliation)
	 sol->upLevelFunction=function->getUPLevelFunction(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP());
          
          return 1;
    }

#include <iostream>
#include <iomanip>      // std::setprecision
using namespace std;
     int LemkeLW::getYLemke(double x[], double y[], double *opt){
       		
		int countNegativeVar=0;
		for(int i=0;i<function->getDimensionLW();i++){
		    if(function->bounds[2*i+2*function->getDimensionUP()]<0)countNegativeVar++;
		}
       
		int sizeLine=(function->getNEQConstraintNumberLW()+function->getDimensionLW() + countNegativeVar);
		int sizeCol=2*(function->getNEQConstraintNumberLW()+function->getDimensionLW() + countNegativeVar)+1+1;
       		double matrix[sizeLine*sizeCol];
		
		//cout<<"\n-------------------------------------------------------\n"<<x[0]<<"\t"<<sizeLine<<"\t"<<sizeCol<<"\n";
		
		//Q,A,c+b
		
		//double xtest[]={0.6,0.39,11.98,17.97};
		
		double q[(function->getDimensionLW()+countNegativeVar)*(function->getDimensionLW()+countNegativeVar)];
		double a[(function->getDimensionLW()+countNegativeVar)*(function->getNEQConstraintNumberLW())];
		double cb[function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW()];
		
		//function->getLemkeMatrix(xtest,matrix,matrix+function->getDimensionLW()*function->getDimensionLW(),matrix+(sizeLine-1)*sizeCol);
		
		function->getLemkeMatrix(x,q,a,cb);
		/*
		for(int i=0;i<function->getDimensionLW()+countNegativeVar;i++){
			for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<q[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
		*/
		if(countNegativeVar>0){
		    for(int i=function->getDimensionLW()-1;i>=0;i--){
		          for(int j=function->getDimensionLW()-1;j>=0;j--){
			  q[i*(function->getDimensionLW()+countNegativeVar)+j]=q[i*function->getDimensionLW()+j];
		          }
		    }
		    
		    for(int i=function->getNEQConstraintNumberLW()-1;i>=0;i--){
		          for(int j=function->getDimensionLW()-1;j>=0;j--){
			  a[i*(function->getDimensionLW()+countNegativeVar)+j]=a[i*function->getDimensionLW()+j];
		          }
		    }
		    
		    for(int i=function->getDimensionLW()+function->getNEQConstraintNumberLW()-1;i>=function->getDimensionLW();i--){
		          cb[i+countNegativeVar]=cb[i];
		    }
		}
		/*
		for(int i=0;i<function->getDimensionLW()+countNegativeVar;i++){
			for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<q[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
		*/
		
		/*
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
			for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<a[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
		*/
		
		/*
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
			for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<a[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";
		*/
		
		for(int i=0, extraVar=0;i<function->getDimensionLW() && extraVar<countNegativeVar;i++){
		    if(function->bounds[2*i+2*function->getDimensionUP()]<0){
		        cb[extraVar+function->getDimensionLW()]=-cb[i];
		        
		        for(int j=0;j<function->getDimensionLW();j++){
			q[j*(function->getDimensionLW()+countNegativeVar)+extraVar+function->getDimensionLW()]=-q[j*(function->getDimensionLW()+countNegativeVar)+i];
			q[(extraVar+function->getDimensionLW())*(function->getDimensionLW()+countNegativeVar)+j]=-q[j*(function->getDimensionLW()+countNegativeVar)+i];
		        }
		        
		        
		        int auxPoint=0;
		        for(int j=0;j<function->getDimensionLW();j++){
			if(function->bounds[2*j+2*function->getDimensionUP()]<0){
			    q[(auxPoint+function->getDimensionLW())*(function->getDimensionLW()+countNegativeVar)+extraVar+function->getDimensionLW()]=q[j*(function->getDimensionLW()+countNegativeVar)+i];
			    auxPoint++;
			}
		        }

		        for(int j=0;j<function->getNEQConstraintNumberLW();j++){
			a[j*(function->getDimensionLW()+countNegativeVar)+extraVar+function->getDimensionLW()]=-a[j*(function->getDimensionLW()+countNegativeVar)+i];
		        }
		        
		        extraVar++;
		    }
		}
		/*
		for(int i=0;i<function->getDimensionLW()+countNegativeVar;i++){
		    for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<q[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";*/
		/*
		for(int i=0;i<function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW();i++){
			cout<<cb[i]<<"\n";
		}
		cout<<"\n\n";
		*/
		/*
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
			for(int j=0;j<function->getDimensionLW()+countNegativeVar;j++) cout<<a[i*(function->getDimensionLW()+countNegativeVar)+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";*/
		
		/*
		matrix[4]=0;
		matrix[5]=0;
		matrix[6]=1;
		matrix[7]=-1;
		cb[0]=1;
		
		matrix[14]=0;
		matrix[15]=0;
		matrix[16]=-1;
		matrix[17]=2;
		cb[1]=4;
		
		matrix[24]=-1;
		matrix[25]=1;
		matrix[26]=2;
		matrix[27]=-2;
		cb[2]=-2;
		
		matrix[34]=1;
		matrix[35]=-2;
		matrix[36]=-2;
		matrix[37]=2;
		cb[3]=-4;*/
		
	
		
		int outBase=0;
		for(int i=0;i<sizeLine;i++){
		    if(cb[i]<cb[outBase]) outBase=i;
		}
		
		if(cb[outBase]>=-TOL_NEQ_CONST){
		   // cout<<"Ja positivo\n";
		   // for(int i=0;i<sizeLine;i++)cout<<cb[i]<<"\n";
		    for(int i=0;i<function->getDimensionLW();i++) y[i]=0;
		    *opt=0;
		    //preencher y
		    return 1;
		}
		
		
		
		
		
		/*
		for(int i=0;i<function->getDimensionLW();i++){
			for(int j=function->getDimensionLW()+function->getNEQConstraintNumberLW();j<2*function->getDimensionLW()+function->getNEQConstraintNumberLW();j++) matrix[i*sizeCol+j]=q[i*function->getDimensionLW()+j-function->getDimensionLW()-function->getNEQConstraintNumberLW()];
		  	
		}
		
		for(int i=0;i<function->getDimensionLW();i++){
		      for(int j=2*function->getDimensionLW()+function->getNEQConstraintNumberLW();j<sizeCol-2;j++) matrix[i*sizeCol+j]=a[i*function->getDimensionLW()+j-2*function->getDimensionLW()-function->getNEQConstraintNumberLW()];
		}
		
		for(int i=function->getDimensionLW();i<sizeLine;i++){
		      for(int j=function->getDimensionLW()+function->getNEQConstraintNumberLW();j<2*function->getDimensionLW()+function->getNEQConstraintNumberLW();j++) matrix[i*sizeCol+j]=a[i*function->getDimensionLW()+j-2*function->getDimensionLW()-function->getNEQConstraintNumberLW()];
		}
		*/
		
		
		for(int i=0;i<sizeLine;i++)
		   for(int j=0;j<sizeCol;j++)matrix[i*sizeCol+j]=0;
		
		//GERANDO M
		for(int i=0;i<function->getDimensionLW()+countNegativeVar;i++){  //Atribuindo Q
		      for(int j=function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW();j<2*(function->getDimensionLW()+countNegativeVar)+function->getNEQConstraintNumberLW();j++) matrix[i*sizeCol+j]=q[i*(function->getDimensionLW()+countNegativeVar)+j-function->getDimensionLW()-countNegativeVar-function->getNEQConstraintNumberLW()];
		  	
		}
		
		for(int i=0;i<function->getDimensionLW()+countNegativeVar;i++){  //Atribuindo A^T
		      for(int j=2*(function->getDimensionLW()+countNegativeVar)+function->getNEQConstraintNumberLW();j<sizeCol-2;j++) matrix[i*sizeCol+j]=a[(j-2*(function->getDimensionLW()+countNegativeVar)-function->getNEQConstraintNumberLW())*(function->getDimensionLW()+countNegativeVar)+i];
		}
		
		for(int i=function->getDimensionLW()+countNegativeVar;i<sizeLine;i++){  //Atribuindo -A
		      for(int j=function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW();j<2*(function->getDimensionLW()+countNegativeVar)+function->getNEQConstraintNumberLW();j++) matrix[i*sizeCol+j]=-a[(i-function->getDimensionLW()-countNegativeVar)*(function->getDimensionLW()+countNegativeVar)+j-function->getDimensionLW()-countNegativeVar-function->getNEQConstraintNumberLW()];
		}
		
		
		
		//Tableau usa -M
		for(int i=0;i<sizeLine;i++)
		    for(int j=0;j<sizeCol;j++) matrix[i*sizeCol+j]*=-1;
		
		for(int i=0;i<sizeLine;i++){
		  	for(int j=0;j<sizeLine;j++) matrix[i*sizeCol+j]=(i==j);
		}
		    
		for(int i=0;i<sizeLine;i++){ //CB com sinais corretos
		      matrix[i*sizeCol+sizeCol-1]=cb[i];
		}   
		      

		for(int i=0;i<sizeLine;i++){
		      matrix[i*sizeCol+sizeCol-2]=-1;
		}
		
		
		for(int i=function->getDimensionLW()+countNegativeVar;i<sizeLine;i++){
		      for(int j=2*(function->getDimensionLW()+countNegativeVar)+function->getNEQConstraintNumberLW();j<sizeCol-2;j++) matrix[i*sizeCol+j]=0;
		}
		
		
		/*
		for(int i=0;i<sizeLine;i++){
			for(int j=0;j<sizeCol;j++) cout<<matrix[i*sizeCol+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";*/

		for(int j=0;j<sizeCol;j++) matrix[outBase*sizeCol+j]*=-1;
		
		
		for(int i=0;i<sizeLine;i++){
			if(i!=outBase)
			    for(int j=0;j<sizeCol;j++) matrix[i*sizeCol+j]+=matrix[outBase*sizeCol+j];
		}
		/*
		for(int i=0;i<sizeLine;i++){
			for(int j=0;j<sizeCol;j++) cout<<matrix[i*sizeCol+j]<<"\t";
			cout<<"\n";
		}
		cout<<"\n\n";*/
		
		int loop=1;
		
		
		int baseVector[function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar];
		for(int i=0;i<function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar;i++)baseVector[i]=i;
		
		baseVector[outBase]=2*(function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar);
		
		int inBase=outBase+function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar;
		
		
		while(loop++){
		  
		        outBase=0;
		        for(int i=0;i<sizeLine;i++){
			if(matrix[sizeCol*i+inBase]>TOL_LEQ_CONST && (matrix[sizeCol*outBase+inBase]<0 || matrix[sizeCol*i+sizeCol-1]/matrix[sizeCol*i+inBase]<=matrix[sizeCol*outBase+sizeCol-1]/matrix[sizeCol*outBase+inBase])){
			    if(matrix[sizeCol*outBase+inBase]<0 || matrix[sizeCol*i+sizeCol-1]/matrix[sizeCol*i+inBase]+TOL_LEQ_CONST<matrix[sizeCol*outBase+sizeCol-1]/matrix[sizeCol*outBase+inBase] || baseVector[i]==2*(function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar))
			      outBase=i;
			}
			//cout<<matrix[sizeCol*i+inBase]<<"\t";
		        }
		        //cout<<"\n";
		        //cout<<"mod base "<<inBase<<"\t"<<outBase<<"\t"<<baseVector[outBase]<<"\n";

		        
		        if(matrix[sizeCol*outBase+inBase]<TOL_LEQ_CONST){
			//Ray termination
		          
			*opt=PENALTY_DEF;
			
			//cout<<"Ray termination\n";
			break;
		        }
		        
		        
		        double multp=matrix[sizeCol*outBase+inBase];
		        for(int j=0;j<sizeCol;j++){
			        matrix[outBase*sizeCol+j]/=multp;
		        }

		        for(int i=0;i<sizeLine;i++){
			        if(i!=outBase){
				        multp=matrix[i*sizeCol+inBase];
				        for(int j=0;j<sizeCol;j++){
					        matrix[i*sizeCol+j]-=multp*matrix[outBase*sizeCol+j];
				        }
			        }
		        }
		        
		        
		        //cout<<"mod base "<<inBase<<"\t"<<outBase<<"\t"<<baseVector[outBase]<<"\n";
	        /*
		        for(int i=0;i<sizeLine;i++){
			        for(int j=0;j<sizeCol;j++) cout<<std::setprecision(2)<<matrix[i*sizeCol+j]<<"\t";
			        cout<<"\n";
		        }
		        cout<<"\n\n";
		       */
		        
		        if(baseVector[outBase]==2*(function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar)){	
			baseVector[outBase]=inBase;
			*opt=0;

			break;
		        }
		        
		        int aux=inBase;
		        if(baseVector[outBase]>=function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar)
			inBase=baseVector[outBase]-function->getNEQConstraintNumberLW()-function->getDimensionLW()-countNegativeVar;
		        else
			inBase=baseVector[outBase]+function->getNEQConstraintNumberLW()+function->getDimensionLW()+countNegativeVar;
		        
		        baseVector[outBase]=aux;
		  
	
		 
		  
		  
		        if(loop>80) {
			        //tableau[function->getDimensionLW()*sizeCol+endCol-1]=1e5;
			        //colocar parada
			        *opt=PENALTY_DEF;
			        cout<<"erro loop";
			        break;
		        }
		}
	
	
	//  double yy[function->getDimensionLW()];
	  for(int i=0;i<function->getDimensionLW();i++)y[i]=0;
	  
	  for(int i=0;i<sizeLine;i++){
	//  cout<<baseVector[i]<<"\n";
	        if(baseVector[i]>=function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW() && baseVector[i]<2*function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW())
		y[baseVector[i] - function->getDimensionLW() - countNegativeVar-function->getNEQConstraintNumberLW()]=matrix[i*sizeCol+sizeCol-1];
	        
	        if(baseVector[i]>=2*function->getDimensionLW()+countNegativeVar+function->getNEQConstraintNumberLW() && baseVector[i]<2*(function->getDimensionLW()+countNegativeVar)+function->getNEQConstraintNumberLW())
		y[baseVector[i] - 2*function->getDimensionLW()- countNegativeVar - function->getNEQConstraintNumberLW()]=-matrix[i*sizeCol+sizeCol-1];
	  }
	  /*
	  cout<<"y:|"<<y[0]<<"|\n";
	  for(int i=1;i<function->getDimensionLW();i++)cout<<"  |"<<y[i]<<"|\n";
	*/
	  return 1;
     }


