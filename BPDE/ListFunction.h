#ifndef LISTFUNCTION_INCLUDED
#define LISTFUNCTION_INCLUDED  
#define DEFINEfunctionListSize 3
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


/* Função 1 - solucao (19,14,0,0.333,0)*/

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
    constraintValuesListReturn[0]=1+(dualNeq[0]*(1) + dualNeq[1]*(5) + dualNeq[2]*(-3)-dualNeq[3]);
    return 1;						
}

const double bound1[4]={0,10e5,0,10e5};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função 2 - solucao (16,11,0,0,3,0,0,0,0)*/

inline double func2UP(double x[], double y[]){  //F(x,y)
  return -x[0]-3*y[0];
}

inline double func2LW(double x[], double y[]){  //f(x,y)
  return -x[0]+3*y[0];
}

inline int func2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int func2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]-2*y[0]+10;
  constraintValuesListReturn[1]=x[0]-2*y[0]-6;
  constraintValuesListReturn[2]=2*x[0]-y[0]-21;
  constraintValuesListReturn[3]=x[0]+2*y[0]-38;
  constraintValuesListReturn[4]=-x[0]+2*y[0]-18;
  constraintValuesListReturn[5]=-x[0];
  constraintValuesListReturn[6]=-y[0];
  return 1;
}


inline int func2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=3+(dualNeq[0]*(-2) + dualNeq[1]*(-2) + dualNeq[2]*(-1) + dualNeq[3]*(2) + dualNeq[4]*(2) + 0*dualNeq[5] - dualNeq[6]);
    return 1;						
}

const double bound2[4]={0,10e5,0,10e5};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função 3 - solucao (19,14,0,0.333,0)*/

inline double func3UP(double x[], double y[]){  //F(x,y)
  return 2*x[0]-11*y[0];
}

inline double func3LW(double x[], double y[]){  //f(x,y)
  return x[0]+3*y[0];
}

inline int func3CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func3CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int func3CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func3CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]-2*y[0]-4;
  constraintValuesListReturn[1]=2*x[0]-y[0]-24;
  constraintValuesListReturn[2]=3*x[0]+4*y[0]-96;
  constraintValuesListReturn[3]=x[0]+7*y[0]-126;
  constraintValuesListReturn[4]=-4*x[0]+5*y[0]-65;
  constraintValuesListReturn[5]=-x[0]-4*y[0]+8;
  constraintValuesListReturn[6]=-x[0];
  constraintValuesListReturn[7]=-y[0];
  return 1;
}


inline int func3CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=3+(dualNeq[0]*(-2) + dualNeq[1]*(-1) + dualNeq[2]*(4) + dualNeq[3]*(7) + dualNeq[4]*(5) + dualNeq[5]*(-4) + 0*dualNeq[6] - dualNeq[7]);
    return 1;						
}

const double bound3[4]={0,10e5,0,10e5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


const Function listFunction[DEFINEfunctionListSize]={{1,1,0,1,0,4,bound1,func1UP,func1LW,func1CTREQUP,func1CTRNEQUP, func1CTREQLW,func1CTRNEQLW,func1CTKKT,"func1"},
						     {1,1,0,0,0,7,bound2,func2UP,func2LW,func2CTREQUP,func2CTRNEQUP, func2CTREQLW,func2CTRNEQLW,func2CTKKT,"func2"},
						     {1,1,0,0,0,8,bound3,func3UP,func3LW,func3CTREQUP,func3CTRNEQUP, func3CTREQLW,func3CTRNEQLW,func3CTKKT,"func3"}};



#endif
