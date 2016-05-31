#ifndef LISTFUNCTIONNEWSMD_INCLUDED
#define LISTFUNCTIONNEWSMD_INCLUDED  
#define DEFINEfunctionListSizeNewSMD 1
#define EPSILON 1e-16
#include <float.h>
#include <cmath>
#include <stdlib.h>



typedef struct functionPrototypeNewSMD{
    int dimensionUP, dimensionLW, numEqConstrUP, numNeqConstrUP, numEqConstrLW, numNeqConstrLW;
    const double *boundsVar;
    double (*funcUP)(double x[], double y[]);
    double (*funcLW)(double x[], double y[]);
    int (*funcCTREQUP)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRNEQUP)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTREQLW)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRNEQLW)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRKKT)(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]);
    int (*funcSimplexTableauKKT)(double x[], double y[],double tableau[]);
    int (*funcLemkeMatrix)(double x[], double matrixQ[], double matrixA[], double matrixCB[]);
    char name[15];
} FunctionNewSMD;

typedef struct functionPrototypeNewSMDIndex{
    int (*setFuncSMD)(FunctionNewSMD **ret);
    char name[15];
} FunctionNewSMDIndex;


extern int inputP, inputQ,inputR,inputS;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD1NEWNewMOD*/

inline double funcSMD1NEWMODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=0;
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=0;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=((y[i+inputQ]*y[i+inputQ])+ 2*x[i+inputP]*x[i+inputP]*x[i+inputP])/inputR; // 2*cos(M_PI*(x[i+inputP]/2.0-(1/2.0))); 
  
  
  return F1+F2+F3;
}

inline double funcSMD1NEWMODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;

  for(int i=0;i<inputQ;i++) f2+=0;
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(x[i+inputP]-y[i+inputQ])*(x[i+inputP]-y[i+inputQ]); 
  
  
  return f1+f2+f3;
}

inline int funcSMD1NEWMODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD1NEWMODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-1;  
      constraintValuesListReturn[2*i+1]=x[i]-1;
  }
  for(int i=0;i<inputR;i++){
      constraintValuesListReturn[2*i+2*(inputP+inputR)]=-((1/4.0)+sin(2*y[i+inputQ])-x[i+inputP]);
      constraintValuesListReturn[2*i+1+2*(inputP+inputR)]=(-1/4.0)+sin(2*y[i+inputQ])-x[i+inputP];
  }
    
  return 1;
}

inline int funcSMD1NEWMODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD1NEWMODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-1;
        constraintValuesListReturn[2*i+1]=y[i]-1;
    }
    
    return 1;
}


inline int funcSMD1NEWMODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=0;

        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(x[i-inputQ+inputP]-y[i])*(-1);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD1NEWMODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
	tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=0;
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(x[i-inputQ+inputP]-y[i])*(-1));
        
    }

    return 1;
}


inline int funcSMD1NEWMODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ;i++){
      
          for(int j=0;j<(inputQ+inputR);j++) matrixQ[i*(inputQ+inputR)+j]=0;

  	  matrixQ[i*(inputQ+inputR) + i]=0; 
      }  

      for(int i=inputQ;i<inputQ+inputR;i++){
      
          for(int j=0;j<(inputQ+inputR);j++) matrixQ[i*(inputQ+inputR)+j]=0;
          
          matrixQ[i*(inputQ+inputR) + i]=2;
      }        

      for(int i=0;i<(inputQ+inputR);i++){
        for(int j=0;j<inputQ+inputR;j++){
	    matrixA[2*i*(inputQ+inputR)+j]=0;
	    matrixA[(2*i+1)*(inputQ+inputR)+j]=0;
	}
	matrixA[2*i*(inputQ+inputR) + i]=-1;
	matrixA[(2*i+1)*(inputQ+inputR) + i]=1;
      }
      
      for(int i=0;i<inputQ;i++){
          matrixCB[i]=0; 
      }
    
      for(int i=inputQ;i<inputQ+inputR;i++){
          matrixCB[i]=-2*(x[i-inputQ+inputP]); 
      }
           
      for(int i=0;i<(inputQ);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }
      for(int i=inputQ;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=1; 
          matrixCB[2*i+inputQ+inputR +1]=1; 
      }


      
      return 1;
}


inline int setFuncSMD1NEWMOD(FunctionNewSMD **ret){
    FunctionNewSMD *func=(FunctionNewSMD *)malloc(sizeof(FunctionNewSMD));
  
    func->funcUP=funcSMD1NEWMODUP;
    func->funcLW=funcSMD1NEWMODLW;
    func->funcCTREQUP=funcSMD1NEWMODCTREQUP;
    func->funcCTRNEQUP=funcSMD1NEWMODCTRNEQUP;
    func->funcCTREQLW=funcSMD1NEWMODCTREQLW;
    func->funcCTRNEQLW=funcSMD1NEWMODCTRNEQLW;
    func->funcCTRKKT=funcSMD1NEWMODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR)+2*inputR;
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD1NEWMOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD1NEWMOD[2*i]=-5;
        boundSMD1NEWMOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD1NEWMOD[2*i]=-1;
        boundSMD1NEWMOD[2*i+1]=1;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD1NEWMOD[2*i]=-5;
        boundSMD1NEWMOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD1NEWMOD[2*i]=-1;
        boundSMD1NEWMOD[2*i+1]=1;
    }

    func->boundsVar=boundSMD1NEWMOD;
    func->funcSimplexTableauKKT=funcSMD1NEWMODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD1NEWMODLemkeMatrix;
    
    (*ret)=func;
    
    return 1;
}

/* --------------------------------------------------------------------------------------------------------------------------------------*/

const FunctionNewSMDIndex listFunctionNewSMD[DEFINEfunctionListSizeNewSMD]={{setFuncSMD1NEWMOD,"funcSMDNew1"},
};



#endif
