#ifndef LISTFUNCTIONSMD_INCLUDED
#define LISTFUNCTIONSMD_INCLUDED  
#define DEFINEfunctionListSizeSMD 3
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
    char name[10];
} FunctionSMD;

typedef struct functionPrototypeSMDIndex{
    int (*setFuncSMD)(FunctionSMD **ret);
    char name[10];
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
  for(int i=0;i<inputR;i++) F3+=(x[i+inputP]*x[i+inputP])+(x[i+inputP]-log10(y[i+inputQ]))*(x[i+inputP]-log10(y[i+inputQ]));
  
  
  return F1+F2+F3;
}

inline double funcSMD2LW(double x[], double y[]){  //f(x,y)
  double f1=0;
  
  for(int i=0;i<inputP;i++) f1+=x[i]*x[i];
  
  double f2=0;
  
  for(int i=0;i<inputQ;i++) f2+=y[i]*y[i];
  
  double f3=0;
  for(int i=0;i<inputR;i++) f3+=(x[i+inputP]-log10(y[i+inputQ]))*(x[i+inputP]-log10(y[i+inputQ]));
  
  
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
        
        constraintValuesListReturn[i]=-2*(x[i-inputQ+inputP]-log10(y[i]))/(y[i]);
       
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
        
        tableau[i*(2*(inputQ+inputR)+1) + 2*(inputQ+inputR)]=2*(x[i-inputQ+inputP]-log10(y[i]))/(y[i]);
        
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
          tableau[i]=-(2*(y[i]-y[i-1]));
        else if(i==inputQ)
          tableau[i]=-(2*(y[i+1]-y[i])*(-1));
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
    
    (*ret)=func;
    
    return 1;
}


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
};



#endif
