#ifndef LISTFUNCTION_INCLUDED
#define LISTFUNCTION_INCLUDED  
#define DEFINEfunctionListSize 17
#include <float.h>
#include <cmath>
#include <stdlib.h>

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
    int (*funcSimplexTableauKKT)(double x[], double y[],double tableau[]);
    int (*funcLemkeMatrix)(double x[], double matrixQ[], double matrixA[], double matrixCB[]);
    double optUPLiterature;
    char name[10];
} Function;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A1 - solucao (19,14,0,0.333,0)*/

inline double funcA1UP(double x[], double y[]){  //F(x,y)
  return x[0]-4*y[0];
}

inline double funcA1LW(double x[], double y[]){  //f(x,y)
  return y[0];
}

inline int funcA1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  return 1;
}

inline int funcA1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-2*x[0]+y[0];
  constraintValuesListReturn[1]=2*x[0]+5*y[0]-108;
  constraintValuesListReturn[2]=2*x[0]-3*y[0]+4;
  constraintValuesListReturn[3]=-y[0];
  return 1;
}


inline int funcA1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=1+(dualNeq[0]*(1) + dualNeq[1]*(5) + dualNeq[2]*(-3)-dualNeq[3]);
    return 1;						
}

inline int funcA1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=1;
	tableau[1]=5;
	tableau[2]=-3;
	tableau[3]=-1;
	tableau[4]=-1;
	return 1;
}

inline int funcA1LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							  //(f)  - (-A  0 )(l1)= (-b)  
      
      matrixA[0]=1;
      
      matrixA[1]=5;
      
      matrixA[2]=-3;
      
      matrixA[3]=-1;
      
      matrixCB[0]=1;
      matrixCB[1]=2*x[0];
      matrixCB[2]=108-2*x[0];
      matrixCB[3]=-4-2*x[0];
      matrixCB[4]=0;

      
      return 1;
}

const double boundA1[4]={0,54,0,18};  //Bounds x, y

//const double boundA1[4]={0,1e3,0,1e3};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A2 - solucao (16,11,0,0,3,0,0,0,0)*/

inline double funcA2UP(double x[], double y[]){  //F(x,y)
  return -x[0]-3*y[0];
}

inline double funcA2LW(double x[], double y[]){  //f(x,y)
  return -x[0]+3*y[0];
}

inline int funcA2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcA2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]-2*y[0]+10;
  constraintValuesListReturn[1]=x[0]-2*y[0]-6;
  constraintValuesListReturn[2]=2*x[0]-y[0]-21;
  constraintValuesListReturn[3]=x[0]+2*y[0]-38;
  constraintValuesListReturn[4]=-x[0]+2*y[0]-18;
  constraintValuesListReturn[5]=-x[0];
  constraintValuesListReturn[6]=-y[0];
  return 1;
}


inline int funcA2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=3+(dualNeq[0]*(-2) + dualNeq[1]*(-2) + dualNeq[2]*(-1) + dualNeq[3]*(2) + dualNeq[4]*(2) + 0*dualNeq[5] - dualNeq[6]);
    return 1;						
}

inline int funcA2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-2;
	tableau[1]=-2;
	tableau[2]=-1;
	tableau[3]=2;
	tableau[4]=2;
	tableau[5]=0;
	tableau[6]=-1;
	tableau[7]=-3;
	return 1;
}


//const double boundA2[4]={0,1e3,0,1e3};  //Bounds x, y

const double boundA2[4]={0,16,0,14};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A3 - solucao (17.45,10.91)*/

inline double funcA3UP(double x[], double y[]){  //F(x,y)
  return 2*x[0]-11*y[0];
}

inline double funcA3LW(double x[], double y[]){  //f(x,y)
  return x[0]+3*y[0];
}

inline int funcA3CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA3CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcA3CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA3CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
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


inline int funcA3CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=3+(dualNeq[0]*(-2) + dualNeq[1]*(-1) + dualNeq[2]*(4) + dualNeq[3]*(7) + dualNeq[4]*(5) + dualNeq[5]*(-4) + 0*dualNeq[6] - dualNeq[7]);
    return 1;						
}

inline int funcA3SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-2;
	tableau[1]=-1;
	tableau[2]=4;
	tableau[3]=7;
	tableau[4]=5;
	tableau[5]=-4;
	tableau[6]=0;
	tableau[7]=-1;
	tableau[8]=-3;
	return 1;
}

//const double boundA3[4]={0,1e3,0,1e3};  //Bounds x, y

const double boundA3[4]={0,18,0,18};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A4 - solucao (0,0.9,0,0.6,0.4)*/

inline double funcA4UP(double x[], double y[]){  //F(x,y)
  return -8*x[0]-4*x[1]+4*y[0]-40*y[1]-4*y[2];
}

inline double funcA4LW(double x[], double y[]){  //f(x,y)
  return x[0]+2*x[1]+y[0]+y[1]+2*y[2];
}

inline int funcA4CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA4CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=-x[1];
  return 1;
}

inline int funcA4CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA4CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0]+y[1]+y[2]-1;
  constraintValuesListReturn[1]=2*x[0]-y[0]+2*y[1]-0.5*y[2]-1;
  constraintValuesListReturn[2]=2*x[1]+2*y[0]-y[1]-0.5*y[2]-1;
  constraintValuesListReturn[3]=-y[0];
  constraintValuesListReturn[4]=-y[1];
  constraintValuesListReturn[5]=-y[2];
  return 1;
}


inline int funcA4CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=1+(dualNeq[0]*(-1) + dualNeq[1]*(-1) + dualNeq[2]*(2) + dualNeq[3]*(-1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(2) + dualNeq[2]*(-1) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[2]=2+(dualNeq[0]*(1) + dualNeq[1]*(-0.5) + dualNeq[2]*(-0.5) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int funcA4SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-1;
	tableau[1]=-1;
	tableau[2]=2;
	tableau[3]=-1;
	tableau[4]=0;
	tableau[5]=0;
	tableau[6]=-1;

	tableau[7]=1;
	tableau[8]=2;
	tableau[9]=-1;
	tableau[10]=0;
	tableau[11]=-1;
	tableau[12]=0;
	tableau[13]=-1;

	tableau[14]=1;
	tableau[15]=-0.5;
	tableau[16]=-0.5;
	tableau[17]=0;
	tableau[18]=0;
	tableau[19]=-1;
	tableau[20]=-2;
	return 1;
}



inline int funcA4LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      matrixQ[4]=0;
      matrixQ[5]=0;
      matrixQ[6]=0;
      matrixQ[7]=0;
      matrixQ[8]=0;
      
      matrixA[0]=-1;
      matrixA[1]=1;
      matrixA[2]=1;
      
      matrixA[3]=-1;
      matrixA[4]=2;
      matrixA[5]=-0.5;
      
      matrixA[6]=2;
      matrixA[7]=-1;
      matrixA[8]=-0.5;
      
      matrixA[9]=-1;
      matrixA[10]=0;
      matrixA[11]=0;
      
      matrixA[12]=0;
      matrixA[13]=-1;
      matrixA[14]=0;
      
      matrixA[15]=0;
      matrixA[16]=0;
      matrixA[17]=-1;
      
      
      matrixCB[0]=1;
      matrixCB[1]=1;
      matrixCB[2]=2;
      matrixCB[3]=1;
      matrixCB[4]=-2*x[0]+1;
      matrixCB[5]=-2*x[1]+1;
      matrixCB[6]=0;
      matrixCB[7]=0;
      matrixCB[8]=0;

      
      return 1;
}


//const double boundA4[10]={0,1e3,0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y
const double boundA4[10]={0,0.75,0,1.5,0,1e3,0,1.5,0,2};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A5 - solucao (1,9,0)*/

inline double funcA5UP(double x[], double y[]){  //F(x,y)
  return -x[0]-2*y[0]-3*y[1];
}

inline double funcA5LW(double x[], double y[]){  //f(x,y)
  return -y[0]+y[1];
}

inline int funcA5CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA5CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-8;
  return 1;
}

inline int funcA5CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA5CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+y[0]+y[1]-10;
  constraintValuesListReturn[1]=-y[0];
  constraintValuesListReturn[2]=y[0]-9;
  constraintValuesListReturn[3]=-y[1];
  constraintValuesListReturn[4]=y[1]-7;
  return 1;
}


inline int funcA5CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-1+(dualNeq[0]*(1) + dualNeq[1]*(-1) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(-1) + dualNeq[4]*(1));
    return 1;						
}

inline int funcA5SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=1;
	tableau[1]=-1;
	tableau[2]=1;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=1;

	tableau[6]=1;
	tableau[7]=0;
	tableau[8]=0;
	tableau[9]=-1;
	tableau[10]=1;
	tableau[11]=-1;

	return 1;
}

const double boundA5[6]={0,8,0,9,0,7};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A6 - solucao (2,0,1.5,0)*/

inline double funcA6UP(double x[], double y[]){  //F(x,y)
  return -2*x[0]+x[1]+0.5*y[0];
}

inline double funcA6LW(double x[], double y[]){  //f(x,y)
  return -4*y[0]+y[1];
}

inline int funcA6CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA6CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]-2;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  return 1;
}

inline int funcA6CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA6CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+y[1]-2.5);
  constraintValuesListReturn[1]=-(-x[0]+3*x[1]-y[1]+2);
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcA6CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-4+(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcA6SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=1;
	tableau[1]=0;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=4;

	tableau[5]=-1;
	tableau[6]=1;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-1;

	return 1;
}

//const double boundA6[8]={0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y
const double boundA6[8]={0,2,0,2,0,5.5,0,8};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A7 - solucao (0.5,0.8,0,0.2,0.8)*/

inline double funcA7UP(double x[], double y[]){  //F(x,y)
  return -8*x[0]-4*x[1]+4*y[0]-40*y[1]-4*y[2];
}

inline double funcA7LW(double x[], double y[]){  //f(x,y)
  return 2*y[0]+y[1]+2*y[2];
}

inline int funcA7CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA7CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+2*x[1]-y[2]-1.3;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  return 1;
}

inline int funcA7CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA7CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0]+y[1]+y[2]-1;
  constraintValuesListReturn[1]=4*x[0]-2*y[0]+4*y[1]-y[2]-2;
  constraintValuesListReturn[2]=4*x[1]+4*y[0]-2*y[1]-y[2]-2;
  constraintValuesListReturn[3]=-y[0];
  constraintValuesListReturn[4]=-y[1];
  constraintValuesListReturn[5]=-y[2];
  return 1;
}


inline int funcA7CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=2+(dualNeq[0]*(-1) + dualNeq[1]*(-2) + dualNeq[2]*(4) + dualNeq[3]*(-1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(4) + dualNeq[2]*(-2) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[2]=2+(dualNeq[0]*(1) + dualNeq[1]*(-1) + dualNeq[2]*(-1) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int funcA7SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-1;
	tableau[1]=-2;
	tableau[2]=4;
	tableau[3]=-1;
	tableau[4]=0;
	tableau[5]=0;
	tableau[6]=-2;

	tableau[7]=1;
	tableau[8]=4;
	tableau[9]=-2;
	tableau[10]=0;
	tableau[11]=-1;
	tableau[12]=0;
	tableau[13]=-1;

	tableau[14]=1;
	tableau[15]=-1;
	tableau[16]=-1;
	tableau[17]=0;
	tableau[18]=0;
	tableau[19]=-1;
	tableau[20]=-2;

	return 1;
}


inline int funcA7LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      matrixQ[4]=0;
      matrixQ[5]=0;
      matrixQ[6]=0;
      matrixQ[7]=0;
      matrixQ[8]=0;

      
      matrixA[0]=-1;
      matrixA[1]=1;
      matrixA[2]=1;
      
      matrixA[3]=-2;
      matrixA[4]=4;
      matrixA[5]=-1;
      
      matrixA[6]=4;
      matrixA[7]=-2;
      matrixA[8]=-1;
      
      
      matrixA[9]=-1;
      matrixA[10]=0;
      matrixA[11]=0;
      
      
      matrixA[12]=0;
      matrixA[13]=-1;
      matrixA[14]=0;
      
      
      matrixA[15]=0;
      matrixA[16]=0;
      matrixA[17]=-1;
      
      matrixCB[0]=2;
      matrixCB[1]=1;
      matrixCB[2]=2;
      matrixCB[3]=1;
      matrixCB[4]=-4*x[0]+2;
      matrixCB[5]=-4*x[1]+2;
      matrixCB[6]=0;
      matrixCB[7]=0;
      matrixCB[8]=0;

      
      return 1;
}



//const double boundA7[10]={0,1e3,0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y

const double boundA7[10]={0,1.5,0,1,0,1.5,0,1.5,0,2};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A8 - solucao (1.55,0.78,0.16,2.21,1.89,0)*/

inline double funcA8UP(double x[], double y[]){  //F(x,y)
  return -4*x[0]+8*x[1]+x[2]-x[3]+9*y[0]-9*y[1];
}

inline double funcA8LW(double x[], double y[]){  //f(x,y)
  return -9*y[0]+9*y[1];
}

inline int funcA8CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA8CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-9*x[0]+3*x[1]-8*x[2]+3*x[3]+3*y[0]-1;
  constraintValuesListReturn[1]=4*x[0]-10*x[1]+3*x[2]+5*x[3]+8*y[0]+8*y[1]-25;
  constraintValuesListReturn[2]=4*x[0]-2*x[1]-2*x[2]+10*x[3]-5*y[0]+8*y[1]-21;
  constraintValuesListReturn[3]=9*x[0]-9*x[1]+4*x[2]-3*x[3]-y[0]-9*y[1]+1;
  constraintValuesListReturn[4]=-2*x[0]-2*x[1]+8*x[2]-5*x[3]+5*y[0]+8*y[1]-20;
  constraintValuesListReturn[5]=7*x[0]+2*x[1]-5*x[2]+4*x[3]-5*y[0]-11;
  constraintValuesListReturn[6]=-x[0];
  constraintValuesListReturn[7]=-x[1];
  constraintValuesListReturn[8]=-x[2];
  constraintValuesListReturn[9]=-x[3];

  return 1;
}

inline int funcA8CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA8CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-6*x[0]+x[1]+x[2]-3*x[3]-9*y[0]-7*y[1]+15;
  constraintValuesListReturn[1]=4*x[1]+5*x[2]+10*x[3]-26;
  constraintValuesListReturn[2]=-9*x[0]+9*x[1]-9*x[2]+5*x[3]-5*y[0]-4*y[1]+5;
  constraintValuesListReturn[3]=5*x[0]+3*x[1]+x[2]+9*x[3]+y[0]+5*y[1]-32;
  constraintValuesListReturn[4]=-y[0];
  constraintValuesListReturn[5]=-y[1];
  return 1;
}


inline int funcA8CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-9+(dualNeq[0]*(-9) + dualNeq[1]*(0) + dualNeq[2]*(-5) + dualNeq[3]*(1) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=9+(dualNeq[0]*(-7) + dualNeq[1]*(0) + dualNeq[2]*(-4) + dualNeq[3]*(5) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int funcA8SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-9;
	tableau[1]=0;
	tableau[2]=-5;
	tableau[3]=1;
	tableau[4]=-1;
	tableau[5]=0;
	tableau[6]=9;

	tableau[7]=-7;
	tableau[8]=0;
	tableau[9]=-4;
	tableau[10]=5;
	tableau[11]=0;
	tableau[12]=-1;
	tableau[13]=-9;

	return 1;
}


inline int funcA8LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      
      matrixA[0]=-9;
      matrixA[1]=-7;
      
      matrixA[2]=0;
      matrixA[3]=0;
      
      matrixA[4]=-5;
      matrixA[5]=-4;
      
      matrixA[6]=1;
      matrixA[7]=5;
      
      matrixA[8]=-1;
      matrixA[9]=-0;
      
      matrixA[10]=0;
      matrixA[11]=-1;
      
   
      matrixCB[0]=-9;
      matrixCB[1]=9;
      matrixCB[2]=6*x[0]-x[1]-x[2]+3*x[3]-15;
      matrixCB[3]=-(4*x[1]+5*x[2]+10*x[3]-26);
      matrixCB[4]=9*x[0]-9*x[1]+9*x[2]-5*x[3]-5;
      matrixCB[5]=-5*x[0]-3*x[1]-x[2]-9*x[3]+32;
      matrixCB[6]=0;
      matrixCB[7]=0;

      
      return 1;
}


const double boundA8[12]={0,5,0,6,0,4,0,3,0,7,0,4};  //Bounds x, y

//const double boundA8[12]={0,1e3,0,1e3,0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

/* Função A9 - solucao (0,2.44,10,0,10,8.74,5.25,10,0,10,3.73,10,10,10,0,0)*/

inline double funcA9UP(double x[], double y[]){  //F(x,y)
  return 12*x[0]-x[1]-12*x[2]+13*x[3]+2*x[5]-5*x[7]+6*x[8]-11*x[9]-5*y[0]-6*y[1]-4*y[2]-7*y[3];
}

inline double funcA9LW(double x[], double y[]){  //f(x,y)
  return 3*y[0]-2*y[1]-3*y[2]-3*y[3]+y[4]+6*y[5];
}

inline int funcA9CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcA9CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  double A[2][10]={{-2,-3,14,-2,-9,2,1,-4,0,2},{1,-7,13,0,-15,2,-8,-4,4,-7}};
  double B[2][6]={{-3,9,-2,-8,1,-8},{-6,-2,6,2,8,-4}};
  double r1[2]={30,-134};

  for(int i=0;i<2;i++){
	constraintValuesListReturn[i]=-r1[i];
	for(int j=0;j<10;j++){
		constraintValuesListReturn[i]+=A[i][j]*x[j];
	}
  }

  for(int i=0;i<2;i++){
	for(int j=0;j<6;j++){
		constraintValuesListReturn[i]+=B[i][j]*y[j];
	}
  }

  for(int i=0;i<10;i++){
	constraintValuesListReturn[2*i+2]=-x[i];
	constraintValuesListReturn[2*i+3]=x[i]-10;
  }

  return 1;
}

inline int funcA9CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcA9CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  double C[7][10]={{-5,7,4,-2,3,-9,9,-1,-3,11},{6,-5,-3,-2,8,5,8,-3,7,3},{-6,-4,2,0,-2,3,-3,2,2,4},{5,6,0,-4,3,-8,1,0,2,-3},{11,-11,4,5,-10,-6,14,-7,-11,-3},{9,-12,-4,-10,2,8,5,-11,-4,1},{7,-2,-6,0,-11,1,-2,-2,-1,-2}};
  double D[7][6]={{10,-9,-6,4,6,-3},{-5,-7,1,1,-6,4},{10,5,6,-4,3,-1},{-4,-3,-4,-4,1,1},{-10,-7,7,7,2,7},{2,-5,10,1,4,5},{-5,-5,-6,-5,1,-12}};

  double r2[7]={83,92,168,-96,-133,89,-192};

 for(int i=0;i<7;i++){
	constraintValuesListReturn[i]=-r2[i];
	for(int j=0;j<10;j++){
		constraintValuesListReturn[i]+=C[i][j]*x[j];
	}
  }

  for(int i=0;i<7;i++){
	for(int j=0;j<6;j++){
		constraintValuesListReturn[i]+=D[i][j]*y[j];
	}
  }

  for(int i=0;i<6;i++){
	constraintValuesListReturn[2*i+7]=-y[i];
	constraintValuesListReturn[2*i+8]=y[i]-10;
  }

  return 1;
}


inline int funcA9CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double D[7][6]={{10,-9,-6,4,6,-3},{-5,-7,1,1,-6,4},{10,5,6,-4,3,-1},{-4,-3,-4,-4,1,1},{-10,-7,7,7,2,7},{2,-5,10,1,4,5},{-5,-5,-6,-5,1,-12}};
    double gradf[6]={3,-2,-3,-3,1,6};

    for(int i=0;i<6;i++){
	constraintValuesListReturn[i]=gradf[i];
	for(int j=0;j<7;j++){
		constraintValuesListReturn[i]+=dualNeq[j]*D[j][i];
	}

	constraintValuesListReturn[i]-=dualNeq[2*i+7];
	constraintValuesListReturn[i]+=dualNeq[2*i+7+1];
    }
 
    return 1;						
}

inline int funcA9SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    double gradf[6]={3,-2,-3,-3,1,6};
    double D[7][6]={{10,-9,-6,4,6,-3},{-5,-7,1,1,-6,4},{10,5,6,-4,3,-1},{-4,-3,-4,-4,1,1},{-10,-7,7,7,2,7},{2,-5,10,1,4,5},{-5,-5,-6,-5,1,-12}};
    for(int i=0;i<6;i++){
	for(int j=0;j<7;j++){
		tableau[i*(19+1)+j]=D[j][i];
	}

	for(int j=7;j<19;j++){
		tableau[i*(19+1)+j]=0;
	}
	
	tableau[i*(19+1) + 2*i + 7]=-1;
	tableau[i*(19+1) + 2*i + 7 + 1]=1;
	tableau[i*(19+1)+19]=-gradf[i];
    }
	

    return 1;
}



inline int funcA9LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      								  //(f)  - (-A  0 )(l1)= (-b)  
      for(int i=0;i<36;i++)
          matrixQ[i]=0;
      
      double D[7][6]={{10,-9,-6,4,6,-3},{-5,-7,1,1,-6,4},{10,5,6,-4,3,-1},{-4,-3,-4,-4,1,1},{-10,-7,7,7,2,7},{2,-5,10,1,4,5},{-5,-5,-6,-5,1,-12}};

      for(int i=0;i<7;i++)
          for(int j=0;j<6;j++)
	  matrixA[i*6+j]=D[i][j];
      
      for(int i=7;i<19;i+=2){
          for(int j=0;j<12;j++)
	  matrixA[i*6+j]=0;
          
          matrixA[i*6 + (i-7)/2]=-1;
          matrixA[(i+1)*6 + (i-7)/2]=1;
      }
	      
      double gradf[6]={3,-2,-3,-3,1,6};
      
      for(int i=0;i<6;i++)   matrixCB[i]=gradf[i];
      
      
      double C[7][10]={{-5,7,4,-2,3,-9,9,-1,-3,11},{6,-5,-3,-2,8,5,8,-3,7,3},{-6,-4,2,0,-2,3,-3,2,2,4},{5,6,0,-4,3,-8,1,0,2,-3},{11,-11,4,5,-10,-6,14,-7,-11,-3},{9,-12,-4,-10,2,8,5,-11,-4,1},{7,-2,-6,0,-11,1,-2,-2,-1,-2}};
      double r2[7]={83,92,168,-96,-133,89,-192};

      for(int i=0;i<7;i++){
          matrixCB[i+6]=r2[i];
          for(int j=0;j<10;j++)
	  matrixCB[i+6]-=C[i][j]*x[j];
      }
      
      for(int i=13;i<25;i+=2){
          matrixCB[i]=0;
          matrixCB[i+1]=10;
      } 
      
      return 1;
}



const double boundA9[32]={0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B1 - solucao (0,30,-10,10)*/

inline double funcB1UP(double x[], double y[]){  //F(x,y)
  return 2*x[0] + 2*x[1] - 3*y[0] - 3*y[1] - 60;
}

inline double funcB1LW(double x[], double y[]){  //f(x,y)
  return (y[0]-x[0]+20)*(y[0]-x[0]+20) + (y[1]-x[1]+20)*(y[1]-x[1]+20);
}

inline int funcB1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]+y[0]-2*y[1]-40;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=x[0]-50;
  constraintValuesListReturn[3]=-x[1];
  constraintValuesListReturn[4]=x[1]-50;

  return 1;
}

inline int funcB1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=2*y[0]-x[0]+10;
  constraintValuesListReturn[1]=2*y[1]-x[1]+10;
  constraintValuesListReturn[2]=-y[0]-10;
  constraintValuesListReturn[3]=y[0]-20;
  constraintValuesListReturn[4]=-y[1]-10;
  constraintValuesListReturn[5]=y[1]-20;

  return 1;
}


inline int funcB1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=2*(y[0]-x[0]+20)+(dualNeq[0]*(2) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=2*(y[1]-x[1]+20)+(dualNeq[0]*(0) + dualNeq[1]*(2) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(1));
  
    return 1;						
}

inline int funcB1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    tableau[0]=2;
    tableau[1]=0;
    tableau[2]=-1;
    tableau[3]=1;
    tableau[4]=0;
    tableau[5]=0;
    tableau[6]=-2*(y[0]-x[0]+20);

    tableau[7]=0;
    tableau[8]=2;
    tableau[9]=0;
    tableau[10]=0;
    tableau[11]=-1;
    tableau[12]=1;
    tableau[13]=-2*(y[1]-x[1]+20);
	

    return 1;
}


inline int funcB1LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=2;

      
      matrixA[0]=2;
      matrixA[1]=0;
      
      matrixA[2]=0;
      matrixA[3]=2;
             
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=1;
      matrixA[7]=0;
      
      matrixA[8]=0;
      matrixA[9]=-1;
      
      matrixA[10]=0;
      matrixA[11]=1;
      
      
      matrixCB[0]=2*(-x[0]+20);
      matrixCB[1]=2*(-x[1]+20);
      matrixCB[2]=x[0]-10;
      matrixCB[3]=x[1]-10;
      matrixCB[4]=10;
      matrixCB[5]=20;
      matrixCB[6]=10;
      matrixCB[7]=20;
      
      
      return 1;
}


const double boundB1[8]={0,50,0,50,-10,20,-10,20};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B2 - solucao (0,2,1.88,0.91)*/

inline double funcB2UP(double x[], double y[]){  //F(x,y)
  return -x[0]*x[0] - 3*x[1] - 4*y[0] + y[1]*y[1];
}

inline double funcB2LW(double x[], double y[]){  //f(x,y)
  return 2*x[0]*x[0]+y[0]*y[0]-5*y[1];
}

inline int funcB2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]*x[0]+2*x[1]-4;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];

  return 1;
}

inline int funcB2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(x[0]*x[0]-2*x[0]+x[1]*x[1]-2*y[0]+y[1]+3);
  constraintValuesListReturn[1]=-(x[1]+3*y[0]-4*y[1]-4);
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];


  return 1;
}


inline int funcB2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=2*y[0]+(dualNeq[0]*(2) + dualNeq[1]*(-3) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=-5+(dualNeq[0]*(-1) + dualNeq[1]*(4) + dualNeq[2]*(0) + dualNeq[3]*(-1));
  
    return 1;						
}

inline int funcB2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    tableau[0]=2;
    tableau[1]=-3;
    tableau[2]=-1;
    tableau[3]=0;
    tableau[4]=-2*(y[0]);


    tableau[5]=-1;
    tableau[6]=4;
    tableau[7]=0;
    tableau[8]=-1;
    tableau[9]=5;

    return 1;
}

//const double boundB2[8]={0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y

const double boundB2[8]={0,2,0,2,0,6,0,4};  //Bounds x, y - Analise Indireta (Não Linear)



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B3 - solucao (7.02, 3.03, 11.98, 17.97, 0.05, 10.00, 29.95, 0)*/

inline double funcB3UP(double x[], double y[]){  //F(x,y)
  return - ( (y[0]+y[2])*(200 - y[0] - y[2]) + (y[1]+y[3])*(160-y[1]-y[3]) ); //MAX=-MIN
}

inline double funcB3LW(double x[], double y[]){  //f(x,y)
  return (y[0]-4)*(y[0]-4)+(y[1]-13)*(y[1]-13)+(y[2]-35)*(y[2]-35)+(y[3]-2)*(y[3]-2);  //MIN= MIN f1 + MIN f2
}

inline int funcB3CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB3CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]+x[2]+x[3]-40;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=x[0]-10;
  constraintValuesListReturn[3]=-x[1];
  constraintValuesListReturn[4]=x[1]-5;
  constraintValuesListReturn[5]=-x[2];
  constraintValuesListReturn[6]=x[2]-15;
  constraintValuesListReturn[7]=-x[3];
  constraintValuesListReturn[8]=x[3]-20;

  return 1;
}

inline int funcB3CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB3CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0         
  constraintValuesListReturn[0]=0.4*y[0] +0.7*y[1] -x[0];
  constraintValuesListReturn[1]=0.6*y[0] +0.3*y[1] -x[1];
  constraintValuesListReturn[2]=0.4*y[2] +0.7*y[3] -x[2];
  constraintValuesListReturn[3]=0.6*y[2] +0.3*y[3] -x[3];
  constraintValuesListReturn[4]=-y[0];
  constraintValuesListReturn[5]=y[0]-20;
  constraintValuesListReturn[6]=-y[1];
  constraintValuesListReturn[7]=y[1]-20;
  constraintValuesListReturn[8]=-y[2];
  constraintValuesListReturn[9]=y[2]-40;
  constraintValuesListReturn[10]=-y[3];
  constraintValuesListReturn[11]=y[3]-40;

  return 1;
}


inline int funcB3CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)

    constraintValuesListReturn[0]=2*(y[0]-4)+(dualNeq[0]*(0.4) + dualNeq[1]*(0.6) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(1) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[1]=2*(y[1]-13)+(dualNeq[0]*(0.7) + dualNeq[1]*(0.3) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(1) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[2]=2*(y[2]-35)+(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(0.4) + dualNeq[3]*(0.6) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(-1) + dualNeq[9]*(1) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[3]=2*(y[3]-2)+(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(0.7) + dualNeq[3]*(0.3) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(-1) + dualNeq[11]*(1));
  
    return 1;						
}

inline int funcB3SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
    tableau[0]=0.4;
    tableau[1]=0.6;
    tableau[2]=0;
    tableau[3]=0;
    tableau[4]=-1;
    tableau[5]=1;
    tableau[6]=0;
    tableau[7]=0;
    tableau[8]=0;
    tableau[9]=0;
    tableau[10]=0;
    tableau[11]=0;
    tableau[12]=-2*(y[0]-4);

    tableau[13]=0.7;
    tableau[14]=0.3;
    tableau[15]=0;
    tableau[16]=0;
    tableau[17]=0;
    tableau[18]=0;
    tableau[19]=-1;
    tableau[20]=1;
    tableau[21]=0;
    tableau[22]=0;
    tableau[23]=0;
    tableau[24]=0;
    tableau[25]=-2*(y[1]-13);


    tableau[26]=0;
    tableau[27]=0;
    tableau[28]=0.4;
    tableau[29]=0.6;
    tableau[30]=0;
    tableau[31]=0;
    tableau[32]=0;
    tableau[33]=0;
    tableau[34]=-1;
    tableau[35]=1;
    tableau[36]=0;
    tableau[37]=0;
    tableau[38]=-2*(y[2]-35);

    tableau[39]=0;
    tableau[40]=0;
    tableau[41]=0.7;
    tableau[42]=0.3;
    tableau[43]=0;
    tableau[44]=0;
    tableau[45]=0;
    tableau[46]=0;
    tableau[47]=0;
    tableau[48]=0;
    tableau[49]=-1;
    tableau[50]=1;
    tableau[51]=-2*(y[3]-2);

    return 1;
}


inline int funcB3LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      
      matrixQ[4]=0;
      matrixQ[5]=2;
      matrixQ[6]=0;
      matrixQ[7]=0;
      
      matrixQ[8]=0;
      matrixQ[9]=0;
      matrixQ[10]=2;
      matrixQ[11]=0;
      
      matrixQ[12]=0;
      matrixQ[13]=0;
      matrixQ[14]=0;
      matrixQ[15]=2;
      
      
      matrixA[0]=0.4;
      matrixA[1]=0.7;
      matrixA[2]=0;
      matrixA[3]=0;
      
      matrixA[4]=0.6;
      matrixA[5]=0.3;
      matrixA[6]=0;
      matrixA[7]=0;
      
      matrixA[8]=0;
      matrixA[9]=0;
      matrixA[10]=0.4;
      matrixA[11]=0.7;
      
      matrixA[12]=0;
      matrixA[13]=0;
      matrixA[14]=0.6;
      matrixA[15]=0.3;
      
      matrixA[16]=-1;
      matrixA[17]=0;
      matrixA[18]=0;
      matrixA[19]=0;
      
      matrixA[20]=1;
      matrixA[21]=0;
      matrixA[22]=0;
      matrixA[23]=0;
      
      matrixA[24]=0;
      matrixA[25]=-1;
      matrixA[26]=0;
      matrixA[27]=0;
      
      matrixA[28]=0;
      matrixA[29]=1;
      matrixA[30]=0;
      matrixA[31]=0;
      
      matrixA[32]=0;
      matrixA[33]=0;
      matrixA[34]=-1;
      matrixA[35]=0;
      
      matrixA[36]=0;
      matrixA[37]=0;
      matrixA[38]=1;
      matrixA[39]=0;
      
      matrixA[40]=0;
      matrixA[41]=0;
      matrixA[42]=0;
      matrixA[43]=-1;
      
      matrixA[44]=0;
      matrixA[45]=0;
      matrixA[46]=0;
      matrixA[47]=1;
      
      
      
      matrixCB[0]=-8;
      matrixCB[1]=-26;
      matrixCB[2]=-70;
      matrixCB[3]=-4;
      matrixCB[4]=x[0];
      matrixCB[5]=x[1];
      matrixCB[6]=x[2];
      matrixCB[7]=x[3];
      matrixCB[8]=0;
      matrixCB[9]=20;
      matrixCB[10]=0;
      matrixCB[11]=20;
      matrixCB[12]=0;
      matrixCB[13]=40;
      matrixCB[14]=0;
      matrixCB[15]=40;

      
      return 1;
}


const double boundB3[16]={0,10,0,5,0,15,0,20,0,20,0,20,0,40,0,40};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B4 - solucao (1.95, 8.05, 0.00, 0.97, 0.97, 1.31, 6.74, 0, 0)*/

inline double funcB4UP(double x[], double y[]){  //F(x,y)
  return -(y[0] * y[1] *sin(x[0]) + y[2] * y[3] * sin(x[1]) + y[4] * y[5] * sin(x[2])); //MIN=-MAX
}

inline double funcB4LW(double x[], double y[]){  //f(x,y)
  return - (y[0]*sin(y[1])+y[1]*sin(y[0]))-(y[2]*sin(y[3])+y[3]*sin(y[2])) -(y[4]*sin(y[5])+y[5]*sin(y[4]));  //MIN= -MAX f1 - MAX f2 - MAX f3
}

inline int funcB4CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB4CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0] + x[1] + x[2] -10;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  constraintValuesListReturn[3]=-x[2];


  return 1;
}

inline int funcB4CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB4CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0       
  constraintValuesListReturn[0]=y[0]+y[1]-x[0];
  constraintValuesListReturn[1]=y[2]+y[3]-x[1];
  constraintValuesListReturn[2]=y[4]+y[5]-x[2];
  constraintValuesListReturn[3]=- y[0];
  constraintValuesListReturn[4]=- y[1];
  constraintValuesListReturn[5]=- y[2];
  constraintValuesListReturn[6]=- y[3];
  constraintValuesListReturn[7]=- y[4];
  constraintValuesListReturn[8]=- y[5];

  return 1;
}


inline int funcB4CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  MAX f(x,y)

    constraintValuesListReturn[0]=(sin(y[1])+y[1]*cos(y[0]))-(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(-1) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0));

    constraintValuesListReturn[1]=(sin(y[0])+y[0]*cos(y[1]))-(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0));

    constraintValuesListReturn[2]=(sin(y[3])+y[3]*cos(y[2]))-(dualNeq[0]*(0) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0));

    constraintValuesListReturn[3]=(sin(y[2])+y[2]*cos(y[3]))-(dualNeq[0]*(0) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(0) + dualNeq[8]*(0));

    constraintValuesListReturn[4]=(sin(y[5])+y[5]*cos(y[4]))-(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(-1) + dualNeq[8]*(0));

    constraintValuesListReturn[5]=(sin(y[4])+y[4]*cos(y[5]))-(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(-1));
  
    return 1;						
}

inline int funcB4SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX f(x,y)
    tableau[0]=1;
    tableau[1]=0;
    tableau[2]=0;
    tableau[3]=-1;
    tableau[4]=0;
    tableau[5]=0;
    tableau[6]=0;
    tableau[7]=0;
    tableau[8]=0;
    tableau[9]=(sin(y[1])+y[1]*cos(y[0]));

    tableau[10]=1;
    tableau[11]=0;
    tableau[12]=0;
    tableau[13]=0;
    tableau[14]=-1;
    tableau[15]=0;
    tableau[16]=0;
    tableau[17]=0;
    tableau[18]=0;
    tableau[19]=(sin(y[0])+y[0]*cos(y[1]));

    tableau[20]=0;
    tableau[21]=1;
    tableau[22]=0;
    tableau[23]=0;
    tableau[24]=0;
    tableau[25]=-1;
    tableau[26]=0;
    tableau[27]=0;
    tableau[28]=0;
    tableau[29]=(sin(y[3])+y[3]*cos(y[2]));

    tableau[30]=0;
    tableau[31]=1;
    tableau[32]=0;
    tableau[33]=0;
    tableau[34]=0;
    tableau[35]=0;
    tableau[36]=-1;
    tableau[37]=0;
    tableau[38]=0;
    tableau[39]=(sin(y[2])+y[2]*cos(y[3]));

    tableau[40]=0;
    tableau[41]=0;
    tableau[42]=1;
    tableau[43]=0;
    tableau[44]=0;
    tableau[45]=0;
    tableau[46]=0;
    tableau[47]=-1;
    tableau[48]=0;
    tableau[49]=(sin(y[5])+y[5]*cos(y[4]));

    tableau[50]=0;
    tableau[51]=0;
    tableau[53]=1;
    tableau[54]=0;
    tableau[55]=0;
    tableau[56]=0;
    tableau[57]=0;
    tableau[58]=-1;
    tableau[59]=(sin(y[4])+y[4]*cos(y[5]));

    return 1;
}

const double boundB4[18]={0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10};  //Bounds x, y

//const double boundB4[18]={0,1e3,0,1e3,0,1e3,0,1e3,0,1e3,0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B5 - solucao (7.07,7.07,7.07,7.07)*/

inline double funcB5UP(double x[], double y[]){  //F(x,y)
  return - ( (x[0] +y[0] )*(x[1] +y[1] )/(1+x[0] *y[0] +x[1] *y[1] )) ; //MIN=-MAX
}

inline double funcB5LW(double x[], double y[]){  //f(x,y)
  return - ( - (x[0] +y[0] )*(x[1] +y[1] )/(1+x[0] *y[0] +x[1] *y[1] ));  //MIN= - MAX (-UP)
}

inline int funcB5CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB5CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0] *x[0] +x[1] *x[1] -100;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];

  return 1;
}

inline int funcB5CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB5CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0    
  constraintValuesListReturn[0]=-y[0];
  constraintValuesListReturn[1]=y[0]-x[0];
  constraintValuesListReturn[2]=-y[1];
  constraintValuesListReturn[3]=y[1]-x[1];

  return 1;
}


inline int funcB5CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  MAX f(x,y)
    constraintValuesListReturn[0]=((x[1]+y[1])*(x[1]*y[1]-x[0]*x[0]+1))/((x[0]*y[0]+x[1]*y[1]+1)*(x[0]*y[0]+x[1]*y[1]+1))-(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=((x[0]+y[0])*(x[0]*y[0]-x[1]*x[1]+1))/((x[0]*y[0]+x[1]*y[1]+1)*(x[0]*y[0]+x[1]*y[1]+1))-(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1));
 
    return 1;						
}

inline int funcB5SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX f(x,y)
    tableau[0]=-1;
    tableau[1]=1;
    tableau[2]=0;
    tableau[3]=0;
    tableau[4]=-((x[1]+y[1])*(x[1]*y[1]-x[0]*x[0]+1))/((x[0]*y[0]+x[1]*y[1]+1)*(x[0]*y[0]+x[1]*y[1]+1));

    tableau[5]=0;
    tableau[6]=0;
    tableau[7]=-1;
    tableau[8]=1;
    tableau[9]=-((x[0]+y[0])*(x[0]*y[0]-x[1]*x[1]+1))/((x[0]*y[0]+x[1]*y[1]+1)*(x[0]*y[0]+x[1]*y[1]+1));

    return 1;
}

//const double boundB5[8]={0,1e3,0,1e3,0,1e3,0,1e3};  //Bounds x, y

const double boundB5[8]={0,10,0,10,0,10,0,10};  //Bounds x, y  - Analise Indireta (Não Linear)



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B6 - solucao (7.07,7.07,7.07,7.07)  IGUALDADE*/ 



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B7 - solucao (0,30,-10,10)*/

inline double funcB7UP(double x[], double y[]){  //F(x,y)
  return fabs(2*x[0]+2*x[1]-3*y[0]-3*y[1]-60) ; 
}

inline double funcB7LW(double x[], double y[]){  //f(x,y)
  return (y[0]-x[0]+20)*(y[0]-x[0]+20) + (y[1]-x[1]+20)*(y[1]-x[1]+20);
}

inline int funcB7CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB7CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]+y[0]-2*y[1]-40;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=x[0]-50;
  constraintValuesListReturn[3]=-x[1];
  constraintValuesListReturn[4]=x[1]-50;

  return 1;
}

inline int funcB7CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB7CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0   
  constraintValuesListReturn[0]=2*y[0]-x[0]+10;
  constraintValuesListReturn[1]=2*y[1]-x[1]+10;
  constraintValuesListReturn[2]=-y[0]-10;
  constraintValuesListReturn[3]=y[0]-20;
  constraintValuesListReturn[4]=-y[1]-10;
  constraintValuesListReturn[5]=y[1]-20;

  return 1;
}


inline int funcB7CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(y[0]-x[0]+20))+(dualNeq[0]*(2) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=(2*(y[1]-x[1]+20))+(dualNeq[0]*(0) + dualNeq[1]*(2) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(1));
 
    return 1;						
}

inline int funcB7SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
    tableau[0]=2;
    tableau[1]=0;
    tableau[2]=-1;
    tableau[3]=1;
    tableau[4]=0;
    tableau[5]=0;
    tableau[6]=-(2*(y[0]-x[0]+20));
    
    tableau[7]=0;
    tableau[8]=2;
    tableau[9]=0;
    tableau[10]=0;
    tableau[11]=-1;
    tableau[12]=1;
    tableau[13]=-(2*(y[1]-x[1]+20));

    return 1;
}


/*Problema y negativo LEMKE*/
inline int funcB7LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=2;

      
      matrixA[0]=2;
      matrixA[1]=0;
      
      matrixA[2]=0;
      matrixA[3]=2;
       
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=1;
      matrixA[7]=0;
      
      matrixA[8]=0;
      matrixA[9]=-1;
      
      matrixA[10]=0;
      matrixA[11]=1;
      
      
      matrixCB[0]=-2*x[0]+40;
      matrixCB[1]=-2*x[1]+40;
      matrixCB[2]=x[0]-10;
      matrixCB[3]=x[1]-10;
      matrixCB[4]=10;
      matrixCB[5]=20;
      matrixCB[6]=10;
      matrixCB[7]=20;
      
      return 1;
}


const double boundB7[8]={0,50,0,50,-10,20,-10,20};  //Bounds x, y




/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B8 - solucao (20,5,10,5)*/

inline double funcB8UP(double x[], double y[]){  //F(x,y)
  return fabs((x[0]-30)*(x[0]-30)+(x[1]-20)*(x[1]-20)-20*y[0]+20*y[1]-225); 
}

inline double funcB8LW(double x[], double y[]){  //f(x,y)
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]);
}

inline int funcB8CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB8CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(x[0]+2*x[1]-30);
  constraintValuesListReturn[1]=x[0]+x[1]-25;
  constraintValuesListReturn[2]=x[1]-15;

  return 1;
}

inline int funcB8CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB8CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0   
  constraintValuesListReturn[0]=-y[0];
  constraintValuesListReturn[1]=y[0]-10;
  constraintValuesListReturn[2]=-y[1]; 
  constraintValuesListReturn[3]=y[1]-10;
  
  return 1;
}


inline int funcB8CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(-(x[0]-y[0])))+(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(2*(-(x[1]-y[1])))+(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1));
 
    return 1;						
}

inline int funcB8SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
    tableau[0]=-1;
    tableau[1]=1;
    tableau[2]=0;
    tableau[3]=0;
    tableau[4]=-(2*(-(x[0]-y[0])));
    
    tableau[5]=0;
    tableau[6]=0;
    tableau[7]=-1;
    tableau[8]=1;
    tableau[9]=-(2*(-(x[1]-y[1])));
    
    return 1;
}


inline int funcB8LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=2;

      
      matrixA[0]=-1;
      matrixA[1]=0;
      
      matrixA[2]=1;
      matrixA[3]=0;
       
      matrixA[4]=0;
      matrixA[5]=-1;
      
      matrixA[6]=0;
      matrixA[7]=1;
      
      
      matrixCB[0]=-2*x[0];
      matrixCB[1]=-2*x[1];
      matrixCB[2]=0;
      matrixCB[3]=10;
      matrixCB[4]=0;
      matrixCB[5]=10;
      
      return 1;
}



//const double boundB8[8]={-1e3,1e3,-1e3,15,0,10,0,10};  //Bounds x, y

const double boundB8[8]={0,20,5,15,0,10,0,10};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B9 - solucao (1.89,0.89,0)*/

inline double funcB9UP(double x[], double y[]){  //F(x,y)
  return fabs((x[0]-1) *(x[0]-1)+2*y[0]-2*x[0]+1.2097); 
}

inline double funcB9LW(double x[], double y[]){  //f(x,y)
  return (2*y[0]-4)*(2*y[0]-4)+(2*y[1]-1)*(2*y[1]-1) + x[0]*y[0];
}

inline int funcB9CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcB9CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];

  return 1;
}

inline int funcB9CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcB9CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0   
  constraintValuesListReturn[0]=4*x[0]+5*y[0]+4*y[1]-12;
  constraintValuesListReturn[1]=4*y[1]-4*x[0]-5*y[0]+4;
  constraintValuesListReturn[2]=4*x[0]-4*y[0]+5*y[1]-4; 
  constraintValuesListReturn[3]=4*y[0]-4*x[0]+5*y[1]-4;
  constraintValuesListReturn[4]=-y[0]; 
  constraintValuesListReturn[5]=-y[1];
  
  return 1;
}


inline int funcB9CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*((2*y[0]-4)+x[0]))+(dualNeq[0]*(5) + dualNeq[1]*(-5) + dualNeq[2]*(-4) + dualNeq[3]*(4) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=(2*((2*y[1]-1)))+(dualNeq[0]*(4) + dualNeq[1]*(4) + dualNeq[2]*(5) + dualNeq[3]*(5) + dualNeq[4]*(0) + dualNeq[5]*(-1));
 
    return 1;						
}

inline int funcB9SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
    tableau[0]=5;
    tableau[1]=-5;
    tableau[2]=-4;
    tableau[3]=4;
    tableau[4]=-1;
    tableau[5]=0;
    tableau[6]=-(2*((2*y[0]-4)+x[0]));
    
    tableau[7]=4;
    tableau[8]=4;
    tableau[9]=5;
    tableau[10]=5;
    tableau[11]=0;
    tableau[12]=-1;
    tableau[13]=-(2*((2*y[1]-1)));
    
    return 1;
}


inline int funcB9LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=4;							  //(f)  - (-A  0 )(l1)= (-b)  
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=4;

      
      matrixA[0]=5;
      matrixA[1]=4;
      
      matrixA[2]=-5;
      matrixA[3]=4;
       
      matrixA[4]=-4;
      matrixA[5]=5;
      
      matrixA[6]=4;
      matrixA[7]=5;
      
      matrixA[8]=-1;
      matrixA[9]=0;
      
      matrixA[10]=0;
      matrixA[11]=-1;
      
      
      matrixCB[0]=x[0]-8;
      matrixCB[1]=-2;
      matrixCB[2]=-4*x[0]+12;
      matrixCB[3]=4*x[0]-4;
      matrixCB[4]=-4*x[0]+4;
      matrixCB[5]=4*x[0]-4;
      matrixCB[6]=0;
      matrixCB[7]=0;
      
      return 1;
}


//const double boundB9[8]={0,1e3,0,1e3,0,1e3};  //Bounds x, y

const double boundB9[6]={0,2,0,2,0,1};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


const Function listFunction[DEFINEfunctionListSize]={{1,1,0,1,0,4,boundA1,funcA1UP,funcA1LW,funcA1CTREQUP,funcA1CTRNEQUP, funcA1CTREQLW,funcA1CTRNEQLW,funcA1CTKKT,funcA1SimplexTableauKKT,funcA1LemkeMatrix,-37,"funcA1"},
				    {1,1,0,0,0,7,boundA2,funcA2UP,funcA2LW,funcA2CTREQUP,funcA2CTRNEQUP, funcA2CTREQLW,funcA2CTRNEQLW,funcA2CTKKT,funcA2SimplexTableauKKT,NULL,-49,"funcA2"},
				    {1,1,0,0,0,8,boundA3,funcA3UP,funcA3LW,funcA3CTREQUP,funcA3CTRNEQUP, funcA3CTREQLW,funcA3CTRNEQLW,funcA3CTKKT,funcA3SimplexTableauKKT,NULL,-85.09,"funcA3"},
				    {2,3,0,2,0,6,boundA4,funcA4UP,funcA4LW,funcA4CTREQUP,funcA4CTRNEQUP, funcA4CTREQLW,funcA4CTRNEQLW,funcA4CTKKT,funcA4SimplexTableauKKT,funcA4LemkeMatrix,-29.2,"funcA4"},
				    {1,2,0,2,0,5,boundA5,funcA5UP,funcA5LW,funcA5CTREQUP,funcA5CTRNEQUP, funcA5CTREQLW,funcA5CTRNEQLW,funcA5CTKKT,funcA5SimplexTableauKKT,NULL,-19,"funcA5"},
				    {2,2,0,3,0,4,boundA6,funcA6UP,funcA6LW,funcA6CTREQUP,funcA6CTRNEQUP, funcA6CTREQLW,funcA6CTRNEQLW,funcA6CTKKT,funcA6SimplexTableauKKT,NULL,-3.25,"funcA6"},
				    {2,3,0,3,0,6,boundA7,funcA7UP,funcA7LW,funcA7CTREQUP,funcA7CTRNEQUP, funcA7CTREQLW,funcA7CTRNEQLW,funcA7CTKKT,funcA7SimplexTableauKKT,funcA7LemkeMatrix,-18.4,"funcA7"},
				    {4,2,0,10,0,6,boundA8,funcA8UP,funcA8LW,funcA8CTREQUP,funcA8CTRNEQUP, funcA8CTREQLW,funcA8CTRNEQLW,funcA8CTKKT,funcA8SimplexTableauKKT,funcA8LemkeMatrix,14.99,"funcA8"},
				    {10,6,0,22,0,19,boundA9,funcA9UP,funcA9LW,funcA9CTREQUP,funcA9CTRNEQUP, funcA9CTREQLW,funcA9CTRNEQLW,funcA9CTKKT,funcA9SimplexTableauKKT,funcA9LemkeMatrix,-453.61,"funcA9"},
				    {2,2,0,5,0,6,boundB1,funcB1UP,funcB1LW,funcB1CTREQUP,funcB1CTRNEQUP, funcB1CTREQLW,funcB1CTRNEQLW,funcB1CTKKT,funcB1SimplexTableauKKT,funcB1LemkeMatrix,0,"funcB1"},
				    {2,2,0,3,0,4,boundB2,funcB2UP,funcB2LW,funcB2CTREQUP,funcB2CTRNEQUP, funcB2CTREQLW,funcB2CTRNEQLW,funcB2CTKKT,funcB2SimplexTableauKKT,NULL,-12.68,"funcB2"},
				    {4,4,0,9,0,12,boundB3,funcB3UP,funcB3LW,funcB3CTREQUP,funcB3CTRNEQUP, funcB3CTREQLW,funcB3CTRNEQLW,funcB3CTKKT,funcB3SimplexTableauKKT,funcB3LemkeMatrix,-6600,"funcB3"},
				    {3,6,0,4,0,9,boundB4,funcB4UP,funcB4LW,funcB4CTREQUP,funcB4CTRNEQUP, funcB4CTREQLW,funcB4CTRNEQLW,funcB4CTKKT,funcB4SimplexTableauKKT,NULL,-9.56,"funcB4"},
				    {2,2,0,3,0,4,boundB5,funcB5UP,funcB5LW,funcB5CTREQUP,funcB5CTRNEQUP, funcB5CTREQLW,funcB5CTRNEQLW,funcB5CTKKT,funcB5SimplexTableauKKT,NULL,-1.98,"funcB5"},
				    {2,2,0,5,0,6,boundB7,funcB7UP,funcB7LW,funcB7CTREQUP,funcB7CTRNEQUP, funcB7CTREQLW,funcB7CTRNEQLW,funcB7CTKKT,funcB7SimplexTableauKKT,funcB7LemkeMatrix,0,"funcB7"},
				    {2,2,0,3,0,4,boundB8,funcB8UP,funcB8LW,funcB8CTREQUP,funcB8CTRNEQUP, funcB8CTREQLW,funcB8CTRNEQLW,funcB8CTKKT,funcB8SimplexTableauKKT,funcB8LemkeMatrix,0,"funcB8"},
				    {1,2,0,1,0,6,boundB9,funcB9UP,funcB9LW,funcB9CTREQUP,funcB9CTRNEQUP, funcB9CTREQLW,funcB9CTRNEQLW,funcB9CTKKT,funcB9SimplexTableauKKT,funcB9LemkeMatrix,0,"funcB9"}
};



#endif
