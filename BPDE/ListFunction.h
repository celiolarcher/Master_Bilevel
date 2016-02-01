#ifndef LISTFUNCTION_INCLUDED
#define LISTFUNCTION_INCLUDED  
#define DEFINEfunctionListSize 2
#include <float.h>


typedef struct functionPrototype{
    int dimensionUP, dimensionLW, numEqConstrUP, numNeqConstrUP, numEqConstrLW, numNeqConstrLW;
    const double *boundsVar;
    double (*funcUP)(double x[], double y[]);
    double (*funcLW)(double x[], double y[]);
    int (*funcCTREQUP)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRNEQUP)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTREQLW)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRNEQLW)(double x[], double y[], double constraintValuesListReturn[]);
    int (*funcCTRKKT)(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]);
    char name[10];
} Function;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função 1 - solucao (19,14,0,0.2,0)*/

inline double func1UP(double x[], double y[]){  //F(x,y)
  return x[0]-4*y[0];
}

inline double func1LW(double x[], double y[]){  //f(x,y)
  return y[0];
}

inline int func1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  return 1;
}

inline int func1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-2*x[0]+y[0];
  constraintValuesListReturn[1]=2*x[0]+5*y[0]-108;
  constraintValuesListReturn[2]=2*x[0]-3*y[0]+4;
  constraintValuesListReturn[3]=-y[0];
  return 1;
}


inline int func1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=1-(dualNeq[0]*(1) + dualNeq[1]*(5) + dualNeq[2]*(-3)-dualNeq[3]);
    return 1;						
}

const double bound1[4]={0,DBL_MAX,0,DBL_MAX};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/



const Function listFunction[DEFINEfunctionListSize]={{1,1,0,1,0,4,bound1,func1UP,func1LW,func1CTREQUP,func1CTRNEQUP, func1CTREQLW,func1CTRNEQLW,func1CTKKT,"func1"},
				     {2,0,1,1,1,1,bound1,func1UP,func1UP,func1CTREQUP,func1CTRNEQUP, func1CTREQLW,func1CTRNEQLW,func1CTKKT,"func2"}};



#endif
