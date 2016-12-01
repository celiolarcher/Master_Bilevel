#ifndef LISTFUNCTIONSMD_INCLUDED
#define LISTFUNCTIONSMD_INCLUDED  
#define DEFINEfunctionListSizeSMD 9
#define EPSILON 1e-16
#include <float.h>
#include <cmath>
#include <stdlib.h>


const double EulerConstant = exp(1.0);

typedef struct functionPrototypeSMD{
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
    double optUPLiterature;
    
} FunctionSMD;

typedef struct functionPrototypeSMDIndex{
    int (*setFuncSMD)(FunctionSMD **ret);
    char name[15];
} FunctionSMDIndex;


extern int inputP, inputQ,inputR,inputS;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD1*/

inline double funcSMD1UP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])+(x[i+inputP]-tan(y[i+inputQ]))*(x[i+inputP]-tan(y[i+inputQ]));
  
  
  return F1+F2+F3;
}

inline double funcSMD1LW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(x[i+inputP]-tan(y[i+inputQ]))*(x[i+inputP]-tan(y[i+inputQ]));
  
  
  return f1+f2+f3;
}

inline int funcSMD1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
    
  return 1;
}

inline int funcSMD1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-M_PI/2;
        constraintValuesListReturn[2*i+1]=y[i]-M_PI/2;
    }
    
    return 1;
}


inline int funcSMD1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        double secy=1/cos(y[i]);
        constraintValuesListReturn[i]=-2*secy*secy*(x[i-inputQ+inputP]-tan(y[i]));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        double secy=1/cos(y[i]);
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=2 * secy*secy *(x[i-inputQ+inputP]-tan(y[i]));
        
    }

    return 1;
}


const double boundSMD1[18]={-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-M_PI/2,M_PI/2,-M_PI/2,M_PI/2};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD2*/

inline double funcSMD2UP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  F2=-F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])+(x[i+inputP]-log(y[i+inputQ]))*(x[i+inputP]-log(y[i+inputQ]));
  
  
  return F1+F2+F3;
}

inline double funcSMD2LW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(x[i+inputP]-log(y[i+inputQ]))*(x[i+inputP]-log(y[i+inputQ]));
  
  
  return f1+f2+f3;
}

inline int funcSMD2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-1;
  }
    
  return 1;
}

inline int funcSMD2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i];
        constraintValuesListReturn[2*i+1]=y[i]-EulerConstant;
    }
    
    return 1;
}


inline int funcSMD2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=-2*(x[i-inputQ+inputP]-log(y[i]))/(y[i]);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=2*(x[i-inputQ+inputP]-log(y[i]))/(y[i]);
        
    }

    return 1;
}


const double boundSMD2[18]={-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,0,EulerConstant,0,EulerConstant};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD6*/

inline double funcSMD6UP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  F2=-F2;
  
  for(int i=inputQ;i<inputQ+inputS;i++) F2+=y[i]*y[i];
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=((x[i+inputP]*x[i+inputP])-(x[i+inputP]-y[i+inputQ+inputS])*(x[i+inputP]-y[i+inputQ+inputS]));
  
  
  return F1+F2+F3;
}

inline double funcSMD6LW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
    
  for(int i=inputQ;i<inputQ+inputS-1;i++) f2+=(y[i+1]-y[i])*(y[i+1]-y[i]);

  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=((x[i+inputP]-y[i+inputQ+inputS])*(x[i+inputP]-y[i+inputQ+inputS]));
  
  
  return f1+f2+f3;
}

inline int funcSMD6CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD6CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
    
  return 1;
}

inline int funcSMD6CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD6CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ+inputS;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ+inputS;i<inputQ+inputS+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    return 1;
}


inline int funcSMD6CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
    for(int i=inputQ;i<inputQ+inputS;i++){
        constraintValuesListReturn[i]=0;
        
        if(i<inputQ+inputS-1)
          constraintValuesListReturn[i]+=2*(y[i+1]-y[i])*(-1);
        if(i>inputQ)
          constraintValuesListReturn[i]+=2*(y[i]-y[i-1]);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ+inputS;i<inputQ+inputS+inputR;i++){
        
        constraintValuesListReturn[i]=2*(x[i-inputQ-inputS+inputP]-y[i])*(-1);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD6SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputS+inputR);j++) tableau[i*(2*(inputQ+inputS+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*(inputQ+inputS+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputS;i++){
      
        for(int j=0;j<2*(inputQ+inputS+inputR);j++) tableau[i*(2*(inputQ+inputS+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i + 1]=1;
        
        
        if(i==inputQ+inputS-1)
          tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*(inputQ+inputS+inputR)]=-(2*(y[i]-y[i-1]));
        else if(i==inputQ)
          tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*(inputQ+inputS+inputR)]=-(2*(y[i+1]-y[i])*(-1));
        else
          tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*(inputQ+inputS+inputR)]=-(2*(y[i+1]-y[i])*(-1) + 2*(y[i]-y[i-1]));
        
    }
    
    for(int i=inputQ+inputS;i<inputQ+inputS+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputS+inputR);j++) tableau[i*(2*(inputQ+inputS+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputS+inputR)+1) + 2*(inputQ+inputS+inputR)]=-(2*(x[i-inputQ-inputS+inputP]-y[i])*(-1));
        
    }

    return 1;
}

inline int funcSMD6LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								    //(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ;i++){
      
          for(int j=0;j<(inputQ+inputS+inputR);j++) matrixQ[i*(inputQ+inputS+inputR)+j]=0;
          
          matrixQ[i*(inputQ+inputS+inputR) + i]=2;
      }
    
      for(int i=inputQ;i<inputQ+inputS;i++){
      
          for(int j=0;j<(inputQ+inputS+inputR);j++) matrixQ[i*(inputQ+inputS+inputR)+j]=0;
          
                   
          if(i==inputQ){
	  matrixQ[i*(inputQ+inputS+inputR) + i]=2; 
	  matrixQ[i*(inputQ+inputS+inputR) + i + 1]=-2; 
          }else if(i==inputQ+inputS-1){
	  matrixQ[i*(inputQ+inputS+inputR) + i]=2; 
	  matrixQ[i*(inputQ+inputS+inputR) + i - 1]=-2; 
          }else{
	  matrixQ[i*(inputQ+inputS+inputR) + i]=4; 
	  matrixQ[i*(inputQ+inputS+inputR) + i - 1]=-2; 
	  matrixQ[i*(inputQ+inputS+inputR) + i + 1]=-2; 
          }
          
      }
    
      for(int i=inputQ+inputS;i<inputQ+inputS+inputR;i++){
      
          for(int j=0;j<(inputQ+inputS+inputR);j++) matrixQ[i*(inputQ+inputS+inputR)+j]=0;
          
          matrixQ[i*(inputQ+inputS+inputR) + i]=2; 
      }
      
      for(int i=0;i<(inputQ+inputS+inputR);i++){
            for(int j=0;j<inputQ+inputS+inputR;j++){
	    matrixA[2*i*(inputQ+inputS+inputR)+j]=0;
	    matrixA[(2*i+1)*(inputQ+inputS+inputR)+j]=0;
	}
	matrixA[2*i*(inputQ+inputS+inputR) + i]=-1;
	matrixA[(2*i+1)*(inputQ+inputS+inputR) + i]=1;
      }
      
      for(int i=0;i<inputQ;i++){
          matrixCB[i]=0; 
      }
    
      for(int i=inputQ;i<inputQ+inputS;i++){
          matrixCB[i]=0; 

      }
    
      for(int i=inputQ+inputS;i<inputQ+inputS+inputR;i++){
          matrixCB[i]=-2*x[i-inputQ-inputS+inputP]; 
      }
           
      for(int i=0;i<(inputQ+inputS+inputR);i++){           
          matrixCB[2*i+inputQ+inputS+inputR]=5; 
          matrixCB[2*i+inputQ+inputS+inputR +1]=10; 
      }


      
      return 1;
}

//const double boundSMD6[20]={-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10,-5,10};  //Bounds x, y


inline int setFuncSMD6(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD6UP;
    func->funcLW=funcSMD6LW;
    func->funcCTREQUP=funcSMD6CTREQUP;
    func->funcCTRNEQUP=funcSMD6CTRNEQUP;
    func->funcCTREQLW=funcSMD6CTREQLW;
    func->funcCTRNEQLW=funcSMD6CTRNEQLW;
    func->funcCTRKKT=funcSMD6CTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputS+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputS+inputR);
    
    double *boundSMD6=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD6[2*i]=-5;
        boundSMD6[2*i+1]=10;
    }
    
    func->boundsVar=boundSMD6;
    func->funcSimplexTableauKKT=funcSMD6SimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD6LemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}




/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD1 MOD*/

inline double funcSMD11UP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  F2-=F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])-(x[i+inputP]-log(y[i+inputQ]))*(x[i+inputP]-log(y[i+inputQ]));
  
  
  return F1+F2+F3;
}

inline double funcSMD11LW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(x[i+inputP]-log(y[i+inputQ]))*(x[i+inputP]-log(y[i+inputQ]));//MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD11CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD11CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-1;
      constraintValuesListReturn[2*i+1]=x[i]-1;
  }
    
  for(int i=0;i<inputR;i++){
      constraintValuesListReturn[2*(inputP+inputR) + i]=1.0/sqrt(inputR)+log(y[i+inputQ])-x[i+inputP];
  }
    
  return 1;
}

inline int funcSMD11CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD11CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]+1/EulerConstant; 
        constraintValuesListReturn[2*i+1]=y[i]-EulerConstant; 
    }
    
    constraintValuesListReturn[2*(inputP+inputR)]=1;
    for(int i=0;i<inputR;i++){
        constraintValuesListReturn[2*(inputP+inputR)]-=(x[i+inputP]-log(y[i+inputQ]))*(x[i+inputP]-log(y[i+inputQ]));
    }
    
    return 1;
}

/*
inline int funcSMD11CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(y[i]-tan(x[i-inputQ+inputP]));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD11SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-tan(x[i-inputQ+inputP])));
        
    }

    return 1;
}


inline int funcSMD11LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ+inputR;i++){
      
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
          matrixCB[i]=-2*tan(x[i-inputQ+inputP]); 
      }
           
      for(int i=0;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }


      
      return 1;
}*/


inline int setFuncSMD11(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD11UP;
    func->funcLW=funcSMD11LW;
    func->funcCTREQUP=funcSMD11CTREQUP;
    func->funcCTRNEQUP=funcSMD11CTRNEQUP;
    func->funcCTREQLW=funcSMD11CTREQLW;
    func->funcCTRNEQLW=funcSMD11CTRNEQLW;
    //func->funcCTRKKT=funcSMD11CTKKT;
    func->funcCTRKKT=NULL;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR) + inputR;
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR) + 1;
    
    double *boundSMD11=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD11[2*i]=-5;
        boundSMD11[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD11[2*i]=-1;
        boundSMD11[2*i+1]=1;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD11[2*i]=-5;
        boundSMD11[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD11[2*i]=1/EulerConstant;
        boundSMD11[2*i+1]=EulerConstant;
    }
    
    func->boundsVar=boundSMD11;
    //func->funcSimplexTableauKKT=funcSMD11SimplexTableauKKT;
    func->funcSimplexTableauKKT=NULL;
    //func->funcLemkeMatrix=funcSMD11LemkeMatrix;
    func->funcLemkeMatrix=NULL;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD1 MOD*/

inline double funcSMD1MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])+(y[i+inputQ]-tan(x[i+inputP]))*(y[i+inputQ]-tan(x[i+inputP]));//MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD1MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(y[i+inputQ]-tan(x[i+inputP]))*(y[i+inputQ]-tan(x[i+inputP]));//MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD1MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD1MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]-M_PI/2+EPSILON;//MODIFICADO
      constraintValuesListReturn[2*i+1]=x[i]-M_PI/2+EPSILON;//MODIFICADO
  }
    
  return 1;
}

inline int funcSMD1MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD1MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;  //MODIFICADO
        constraintValuesListReturn[2*i+1]=y[i]-10;  //MODIFICADO
    }
    
    return 1;
}


inline int funcSMD1MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(y[i]-tan(x[i-inputQ+inputP]));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD1MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-tan(x[i-inputQ+inputP])));
        
    }

    return 1;
}


inline int funcSMD1MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ+inputR;i++){
      
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
          matrixCB[i]=-2*tan(x[i-inputQ+inputP]); 
      }
           
      for(int i=0;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }


      
      return 1;
}


inline int setFuncSMD1MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD1MODUP;
    func->funcLW=funcSMD1MODLW;
    func->funcCTREQUP=funcSMD1MODCTREQUP;
    func->funcCTRNEQUP=funcSMD1MODCTRNEQUP;
    func->funcCTREQLW=funcSMD1MODCTREQLW;
    func->funcCTRNEQLW=funcSMD1MODCTRNEQLW;
    func->funcCTRKKT=funcSMD1MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD1MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD1MOD[2*i]=-5;
        boundSMD1MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD1MOD[2*i]=-M_PI/2+EPSILON;
        boundSMD1MOD[2*i+1]=M_PI/2-EPSILON;
    }

    for(int i=inputP+inputR;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD1MOD[2*i]=-5;
        boundSMD1MOD[2*i+1]=10;
    }
    
    func->boundsVar=boundSMD1MOD;
    func->funcSimplexTableauKKT=funcSMD1MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD1MODLemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD2MOD*/

inline double funcSMD2MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  F2=-F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]-1)*(x[i+inputP]-1)-(y[i+inputQ]-log(x[i+inputP]))*(y[i+inputQ]-log(x[i+inputP])); //MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD2MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(y[i+inputQ]-log(x[i+inputP]))*(y[i+inputQ]-log(x[i+inputP]));
  
  
  return f1+f2+f3;
}

inline int funcSMD2MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD2MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]+1.0/EulerConstant;//MODIFICADO
      constraintValuesListReturn[2*i+1]=x[i]-EulerConstant;//MODIFICADO
  }
    
  return 1;
}

inline int funcSMD2MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD2MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;//MODIFICADO
        constraintValuesListReturn[2*i+1]=y[i]-1;//MODIFICADO
    }
    
    return 1;
}


inline int funcSMD2MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(y[i]-log(x[i-inputQ+inputP]));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD2MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-log(x[i-inputQ+inputP])));
        
    }

    return 1;
}


inline int funcSMD2MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ+inputR;i++){
      
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
          matrixCB[i]=-2*log(x[i-inputQ+inputP]); 
      }
           
      for(int i=0;i<(inputQ);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }
      for(int i=inputQ;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=1; 
      }


      
      return 1;
}


inline int setFuncSMD2MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD2MODUP;
    func->funcLW=funcSMD2MODLW;
    func->funcCTREQUP=funcSMD2MODCTREQUP;
    func->funcCTRNEQUP=funcSMD2MODCTRNEQUP;
    func->funcCTREQLW=funcSMD2MODCTREQLW;
    func->funcCTRNEQLW=funcSMD2MODCTRNEQLW;
    func->funcCTRKKT=funcSMD2MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD2MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD2MOD[2*i]=-5;
        boundSMD2MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD2MOD[2*i]=1.0/EulerConstant;
        boundSMD2MOD[2*i+1]=EulerConstant;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD2MOD[2*i]=-5;
        boundSMD2MOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD2MOD[2*i]=-5;
        boundSMD2MOD[2*i+1]=1;
    }

    func->boundsVar=boundSMD2MOD;
    func->funcSimplexTableauKKT=funcSMD2MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD2MODLemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD3MOD*/

inline double funcSMD3MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])+(y[i+inputQ]-tan(sqrt(x[i+inputP])))*(y[i+inputQ]-tan(sqrt(x[i+inputP]))); //MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD3MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=inputQ;
  
  for(int i=0;i<inputQ;i++) f2+=((y[i]*y[i])-1);
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(y[i+inputQ]-tan(sqrt(x[i+inputP])))*(y[i+inputQ]-tan(sqrt(x[i+inputP]))); //MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD3MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD3MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i];//MODIFICADO
      constraintValuesListReturn[2*i+1]=x[i]-M_PI/2+EPSILON;//MODIFICADO
  }
    
  return 1;
}

inline int funcSMD3MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD3MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;//MODIFICADO
        constraintValuesListReturn[2*i+1]=y[i]-10;//MODIFICADO
    }
    
    return 1;
}


inline int funcSMD3MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(y[i]-tan(sqrt(x[i-inputQ+inputP])));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD3MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-tan(sqrt(x[i-inputQ+inputP]))));
        
    }

    return 1;
}


inline int funcSMD3MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ+inputR;i++){
      
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
          matrixCB[i]=-2*tan(sqrt(x[i-inputQ+inputP])); 
      }
           
      for(int i=0;i<(inputQ);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }
      for(int i=inputQ;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }


      
      return 1;
}


inline int setFuncSMD3MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD3MODUP;
    func->funcLW=funcSMD3MODLW;
    func->funcCTREQUP=funcSMD3MODCTREQUP;
    func->funcCTRNEQUP=funcSMD3MODCTRNEQUP;
    func->funcCTREQLW=funcSMD3MODCTREQLW;
    func->funcCTRNEQLW=funcSMD3MODCTRNEQLW;
    func->funcCTRKKT=funcSMD3MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD3MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD3MOD[2*i]=-5;
        boundSMD3MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD3MOD[2*i]=0;
        boundSMD3MOD[2*i+1]=M_PI/2-EPSILON;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD3MOD[2*i]=-5;
        boundSMD3MOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD3MOD[2*i]=-5;
        boundSMD3MOD[2*i+1]=10;
    }

    func->boundsVar=boundSMD3MOD;
    func->funcSimplexTableauKKT=funcSMD3MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD3MODLemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD4MOD*/

inline double funcSMD4MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];

  F2=-F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])-(y[i+inputQ]-log(1+x[i+inputP]))*(y[i+inputQ]-log(1+x[i+inputP])); //MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD4MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=inputQ;
  
  for(int i=0;i<inputQ;i++) f2+=((y[i]*y[i])-1);
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(y[i+inputQ]-log(1+x[i+inputP]))*(y[i+inputQ]-log(1+x[i+inputP])); //MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD4MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD4MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i];//MODIFICADO
      constraintValuesListReturn[2*i+1]=x[i]-EulerConstant;//MODIFICADO
  }
    
  return 1;
}

inline int funcSMD4MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD4MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-1;//MODIFICADO
        constraintValuesListReturn[2*i+1]=y[i]-1;//MODIFICADO
    }
    
    return 1;
}


inline int funcSMD4MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(y[i]-log(1+x[i-inputQ+inputP]));
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD4MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-2*y[i];
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-log(1+x[i-inputQ+inputP])));
        
    }

    return 1;
}


inline int funcSMD4MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ+inputR;i++){
      
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
          matrixCB[i]=-2*log(1+x[i-inputQ+inputP]); 
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


inline int setFuncSMD4MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD4MODUP;
    func->funcLW=funcSMD4MODLW;
    func->funcCTREQUP=funcSMD4MODCTREQUP;
    func->funcCTRNEQUP=funcSMD4MODCTRNEQUP;
    func->funcCTREQLW=funcSMD4MODCTREQLW;
    func->funcCTRNEQLW=funcSMD4MODCTRNEQLW;
    func->funcCTRKKT=funcSMD4MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD4MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD4MOD[2*i]=-5;
        boundSMD4MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD4MOD[2*i]=0;
        boundSMD4MOD[2*i+1]=EulerConstant;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD4MOD[2*i]=-5;
        boundSMD4MOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD4MOD[2*i]=-1;
        boundSMD4MOD[2*i+1]=1;
    }

    func->boundsVar=boundSMD4MOD;
    func->funcSimplexTableauKKT=funcSMD4MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD4MODLemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD5MOD*/

inline double funcSMD5MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ-1;i++) F2+=(y[i]-1)*(y[i]-1);
    
  for(int i=0;i<inputQ-1;i++) F2+=(y[i+1]-y[i])*(y[i+1]-y[i]);  //MODIFICADO
  
  F2=-F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])-(sqrt(x[i+inputP])-y[i+inputQ])*(sqrt(x[i+inputP])-y[i+inputQ]); //MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD5MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;

  for(int i=0;i<inputQ-1;i++) f2+=(y[i]-1)*(y[i]-1);
    
  for(int i=0;i<inputQ-1;i++) f2+=(y[i+1]-y[i])*(y[i+1]-y[i]);  //MODIFICADO
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(sqrt(x[i+inputP])-y[i+inputQ])*(sqrt(x[i+inputP])-y[i+inputQ]); //MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD5MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD5MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i];  //MODIFICADO >=0
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
    
  return 1;
}

inline int funcSMD5MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD5MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    return 1;
}


inline int funcSMD5MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
	if(i==0){
	        constraintValuesListReturn[i]=2*(y[i]-1)+2*(y[i+1]-y[i])*(-1);
        }else if(i==inputQ-1){
	        constraintValuesListReturn[i]=2*(y[i]-y[i-1]);
	}else{
	        constraintValuesListReturn[i]=2*(y[i]-1)+2*(y[i+1]-y[i])*(-1)+2*(y[i]-y[i-1]);
	}

        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(sqrt(x[i-inputQ+inputP])-y[i])*(-1);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD5MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
	if(i==0){
		tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-1)+2*(y[i+1]-y[i])*(-1));
        }else if(i==inputQ-1){
        	tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-y[i-1]));
	}else{
        	tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(y[i]-1)+2*(y[i+1]-y[i])*(-1)+2*(y[i]-y[i-1]));
	}
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(sqrt(x[i-inputQ+inputP])-y[i])*(-1));
        
    }

    return 1;
}


inline int funcSMD5MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ;i++){
      
          for(int j=0;j<(inputQ+inputR);j++) matrixQ[i*(inputQ+inputR)+j]=0;



 	  if(i==0){
	  	  matrixQ[i*(inputQ+inputR) + i]=4; 
	  	  matrixQ[i*(inputQ+inputR) + i + 1]=-2; 
          }else if(i==inputQ-1){
		  matrixQ[i*(inputQ+inputR) + i]=2; 
		  matrixQ[i*(inputQ+inputR) + i - 1]=-2; 
          }else{
		  matrixQ[i*(inputQ+inputR) + i]=6; 
		  matrixQ[i*(inputQ+inputR) + i - 1]=-2; 
		  matrixQ[i*(inputQ+inputR) + i + 1]=-2; 
          }
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
      
      for(int i=0;i<inputQ-1;i++){
          matrixCB[i]=-2; 
      }
      
      if(inputQ>0)
          matrixCB[inputQ-1]=0; 

      for(int i=inputQ;i<inputQ+inputR;i++){
          matrixCB[i]=-2*(sqrt(x[i-inputQ+inputP])); 
      }
           
      for(int i=0;i<(inputQ);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }
      for(int i=inputQ;i<(inputQ+inputR);i++){           
          matrixCB[2*i+inputQ+inputR]=5; 
          matrixCB[2*i+inputQ+inputR +1]=10; 
      }


      
      return 1;
}


inline int setFuncSMD5MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD5MODUP;
    func->funcLW=funcSMD5MODLW;
    func->funcCTREQUP=funcSMD5MODCTREQUP;
    func->funcCTRNEQUP=funcSMD5MODCTRNEQUP;
    func->funcCTREQLW=funcSMD5MODCTREQLW;
    func->funcCTRNEQLW=funcSMD5MODCTRNEQLW;
    func->funcCTRKKT=funcSMD5MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR);
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD5MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD5MOD[2*i]=-5;
        boundSMD5MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD5MOD[2*i]=0;
        boundSMD5MOD[2*i+1]=10;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD5MOD[2*i]=-5;
        boundSMD5MOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD5MOD[2*i]=-5;
        boundSMD5MOD[2*i+1]=10;
    }

    func->boundsVar=boundSMD5MOD;
    func->funcSimplexTableauKKT=funcSMD5MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD5MODLemkeMatrix;
    
    func->optUPLiterature=0;
    
    (*ret)=func;
    
    return 1;
}


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função SMD11MOD*/
/*
inline double funcSMD11MODUP(double x[], double y[]){  //F(x,y)
  double F1=0;
  
  for(int i=0;i<inputP;i++) F1+=x[i]*x[i];
  
  double F2=0;
  
  for(int i=0;i<inputQ;i++) F2+=y[i]*y[i];
  
  F2=-F2;
  
  double F3=0;
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])-(log(x[i+inputP])-y[i+inputQ])*(log(x[i+inputP])-y[i+inputQ]); //MODIFICADO
  
  
  return F1+F2+F3;
}

inline double funcSMD11MODLW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;

  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(log(x[i+inputP])-y[i+inputQ])*(log(x[i+inputP])-y[i+inputQ]); //MODIFICADO
  
  
  return f1+f2+f3;
}

inline int funcSMD11MODCTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcSMD11MODCTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  for(int i=0;i<inputP;i++){
      constraintValuesListReturn[2*i]=-x[i]-5;
      constraintValuesListReturn[2*i+1]=x[i]-10;
  }
  
  for(int i=inputP;i<inputP+inputR;i++){
      constraintValuesListReturn[2*i]=-x[i]+1/EulerConstant;  //MODIFICADO 
      constraintValuesListReturn[2*i+1]=x[i]-EulerConstant;
  }

  for(int i=0;i<inputR;i++){
	constraintValuesListReturn[2*(inputP+inputR)+i]=-y[i+inputQ]+1/(sqrt(inputR))+log(x[inputP+i]);
  }
    
  return 1;
}

inline int funcSMD11MODCTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcSMD11MODCTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[2*i]=-y[i]-5;
        constraintValuesListReturn[2*i+1]=y[i]-10;
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
        constraintValuesListReturn[2*i]=-y[i]-1; //MODIFICADO
        constraintValuesListReturn[2*i+1]=y[i]-1;  //MODIFICADO
    }
    
    return 1;
}


inline int funcSMD11MODCTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)   
    for(int i=0;i<inputQ;i++){
        constraintValuesListReturn[i]=2*y[i];

        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
  
  
    for(int i=inputQ;i<inputQ+inputR;i++){
        
        constraintValuesListReturn[i]=2*(log(x[i-inputQ+inputP])-y[i])*(-1);
       
        constraintValuesListReturn[i]-=dualNeq[2*i];
        constraintValuesListReturn[i]+=dualNeq[2*i+1];
    }
    
    return 1;						
}

inline int funcSMD11MODSimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    for(int i=0;i<inputQ;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
	tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*y[i]);
        
    }
    
    for(int i=inputQ;i<inputQ+inputR;i++){
      
        for(int j=0;j<2*(inputQ+inputR);j++) tableau[i*(2*(inputQ+inputR)+1)+j]=0;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*i]=-1;
        tableau[i*(2*(inputQ+inputR)+1) + 2*i + 1]=1;
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=-(2*(log(x[i-inputQ+inputP])-y[i])*(-1));
        
    }

    return 1;
}


inline int funcSMD11MODLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
								 					//(f)  - (-A  0 )(l1)= (-b)  
     
      for(int i=0;i<inputQ;i++){
      
          for(int j=0;j<(inputQ+inputR);j++) matrixQ[i*(inputQ+inputR)+j]=0;

  	  matrixQ[i*(inputQ+inputR) + i]=2; 
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
          matrixCB[i]=-2*(log(x[i-inputQ+inputP])); 
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


inline int setFuncSMD11MOD(FunctionSMD **ret){
    FunctionSMD *func=(FunctionSMD *)malloc(sizeof(FunctionSMD));
  
    func->funcUP=funcSMD11MODUP;
    func->funcLW=funcSMD11MODLW;
    func->funcCTREQUP=funcSMD11MODCTREQUP;
    func->funcCTRNEQUP=funcSMD11MODCTRNEQUP;
    func->funcCTREQLW=funcSMD11MODCTREQLW;
    func->funcCTRNEQLW=funcSMD11MODCTRNEQLW;
    func->funcCTRKKT=funcSMD11MODCTKKT;
    func->dimensionUP=inputP+inputR;
    func->dimensionLW=inputQ+inputR;
    func->numEqConstrUP=0;
    func->numNeqConstrUP=2*(inputP+inputR) + inputR;
    func->numEqConstrLW=0;
    func->numNeqConstrLW=2*(inputQ+inputR);
    
    double *boundSMD11MOD=new double[2*(func->dimensionUP+func->dimensionLW)];
    
    for(int i=0;i<inputP;i++){
        boundSMD11MOD[2*i]=-5;
        boundSMD11MOD[2*i+1]=10;
    }

    for(int i=inputP;i<inputP+inputR;i++){
        boundSMD11MOD[2*i]=EulerConstant;
        boundSMD11MOD[2*i+1]=1/EulerConstant;
    }

    for(int i=inputP+inputR;i<inputP+inputR+inputQ;i++){
        boundSMD11MOD[2*i]=-5;
        boundSMD11MOD[2*i+1]=10;
    }
    
    for(int i=inputP+inputR+inputQ;i<func->dimensionUP+func->dimensionLW;i++){
        boundSMD11MOD[2*i]=-1;
        boundSMD11MOD[2*i+1]=1;
    }

    func->boundsVar=boundSMD11MOD;
    func->funcSimplexTableauKKT=funcSMD11MODSimplexTableauKKT;
    func->funcLemkeMatrix=funcSMD11MODLemkeMatrix;
    
    (*ret)=func;
    
    return 1;
}
*/
/* --------------------------------------------------------------------------------------------------------------------------------------*/

/*
const FunctionSMD listFunctionSMD[DEFINEfunctionListSizeSMD]={{5,4,0,10,0,8,boundSMD1,funcSMD1UP,funcSMD1LW,funcSMD1CTREQUP,funcSMD1CTRNEQUP, funcSMD1CTREQLW,funcSMD1CTRNEQLW,funcSMD1CTKKT,funcSMD1SimplexTableauKKT,NULL,"funcSMD1"},
					  {5,4,0,10,0,8,boundSMD2,funcSMD2UP,funcSMD2LW,funcSMD2CTREQUP,funcSMD2CTRNEQUP, funcSMD2CTREQLW,funcSMD2CTRNEQLW,funcSMD2CTKKT,funcSMD2SimplexTableauKKT,NULL,"funcSMD2"},
					  {5,5,0,10,0,10,boundSMD6,funcSMD6UP,funcSMD6LW,funcSMD6CTREQUP,funcSMD6CTRNEQUP, funcSMD6CTREQLW,funcSMD6CTRNEQLW,funcSMD6CTKKT,funcSMD6SimplexTableauKKT,funcSMD6LemkeMatrix,"funcSMD6"},
};
*/
const FunctionSMDIndex listFunctionSMD[DEFINEfunctionListSizeSMD]={{NULL,"funcSMD1"},
					       {NULL,"funcSMD2"},
					       {setFuncSMD6,"funcSMD6"},
					       {setFuncSMD11,"funcSMD11"},
					       {setFuncSMD1MOD,"funcSMD1MOD"},
					       {setFuncSMD2MOD,"funcSMD2MOD"},
					       {setFuncSMD3MOD,"funcSMD3MOD"},
					       {setFuncSMD4MOD,"funcSMD4MOD"},
					       {setFuncSMD5MOD,"funcSMD5MOD"},
						//{setFuncSMD11MOD,"funcSMD11MOD"},
};



#endif
