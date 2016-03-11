#ifndef LISTFUNCTION_INCLUDED
#define LISTFUNCTION_INCLUDED  
#define DEFINEfunctionListSize 14
#include <float.h>
#include <cmath>

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
    char name[10];
} Function;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A1 - solucao (19,14,0,0.333,0)*/

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

inline int func1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=1;
	tableau[1]=5;
	tableau[2]=-3;
	tableau[3]=-1;
	tableau[4]=-1;
	return 1;
}

const double bound1[4]={0,10e5,0,10e5};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A2 - solucao (16,11,0,0,3,0,0,0,0)*/

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

inline int func2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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


const double bound2[4]={0,10e5,0,10e5};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A3 - solucao (17.45,10.91)*/

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

inline int func3SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound3[4]={0,10e5,0,10e5};  //Bounds x, y




/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A4 - solucao (0,0.9,0,0.6,0.4)*/

inline double func4UP(double x[], double y[]){  //F(x,y)
  return -8*x[0]-4*x[1]+4*y[0]-40*y[1]-4*y[2];
}

inline double func4LW(double x[], double y[]){  //f(x,y)
  return x[0]+2*x[1]+y[0]+y[1]+2*y[2];
}

inline int func4CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func4CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=-x[1];
  return 1;
}

inline int func4CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func4CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0]+y[1]+y[2]-1;
  constraintValuesListReturn[1]=2*x[0]-y[0]+2*y[1]-0.5*y[2]-1;
  constraintValuesListReturn[2]=2*x[1]+2*y[0]-y[1]-0.5*y[2]-1;
  constraintValuesListReturn[3]=-y[0];
  constraintValuesListReturn[4]=-y[1];
  constraintValuesListReturn[5]=-y[2];
  return 1;
}


inline int func4CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=1+(dualNeq[0]*(-1) + dualNeq[1]*(-1) + dualNeq[2]*(2) + dualNeq[3]*(-1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(2) + dualNeq[2]*(-1) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[2]=2+(dualNeq[0]*(1) + dualNeq[1]*(-0.5) + dualNeq[2]*(-0.5) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int func4SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound4[10]={0,10e5,0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A5 - solucao (1,9,0)*/

inline double func5UP(double x[], double y[]){  //F(x,y)
  return -x[0]-2*y[0]-3*y[1];
}

inline double func5LW(double x[], double y[]){  //f(x,y)
  return -y[0]+y[1];
}

inline int func5CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func5CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-8;
  return 1;
}

inline int func5CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func5CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+y[0]+y[1]-10;
  constraintValuesListReturn[1]=-y[0];
  constraintValuesListReturn[2]=y[0]-9;
  constraintValuesListReturn[3]=-y[1];
  constraintValuesListReturn[4]=y[1]-7;
  return 1;
}


inline int func5CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-1+(dualNeq[0]*(1) + dualNeq[1]*(-1) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(-1) + dualNeq[4]*(1));
    return 1;						
}

inline int func5SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound5[6]={0,8,0,9,0,7};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A6 - solucao (2,0,1.5,0)*/

inline double func6UP(double x[], double y[]){  //F(x,y)
  return -2*x[0]+x[1]+0.5*y[0];
}

inline double func6LW(double x[], double y[]){  //f(x,y)
  return -4*y[0]+y[1];
}

inline int func6CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func6CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]-2;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  return 1;
}

inline int func6CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func6CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+y[1]-2.5);
  constraintValuesListReturn[1]=-(-x[0]+3*x[1]-y[1]+2);
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int func6CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-4+(dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int func6SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound6[8]={0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A7 - solucao (0.5,0.8,0,0.2,0.8)*/

inline double func7UP(double x[], double y[]){  //F(x,y)
  return -8*x[0]-4*x[1]+4*y[0]-40*y[1]-4*y[2];
}

inline double func7LW(double x[], double y[]){  //f(x,y)
  return 2*y[0]+y[1]+y[2];
}

inline int func7CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func7CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+2*x[1]-y[2]-1.3;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  return 1;
}

inline int func7CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func7CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0]+y[1]+y[2]-1;
  constraintValuesListReturn[1]=4*x[0]-2*y[0]+4*y[1]-y[2]-2;
  constraintValuesListReturn[2]=4*x[1]+4*y[0]-2*y[1]-y[2]-2;
  constraintValuesListReturn[3]=-y[0];
  constraintValuesListReturn[4]=-y[1];
  constraintValuesListReturn[5]=-y[2];
  return 1;
}


inline int func7CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=2+(dualNeq[0]*(-1) + dualNeq[1]*(-2) + dualNeq[2]*(4) + dualNeq[3]*(-1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=1+(dualNeq[0]*(1) + dualNeq[1]*(4) + dualNeq[2]*(-2) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[2]=1+(dualNeq[0]*(1) + dualNeq[1]*(-1) + dualNeq[2]*(-1) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int func7SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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
	tableau[20]=-1;

	return 1;
}

const double bound7[10]={0,10e5,0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função A8 - solucao (1.55,0.78,0.16,2.21,1.89,0)*/

inline double func8UP(double x[], double y[]){  //F(x,y)
  return -4*x[0]+8*x[1]+x[2]-x[3]+9*y[0]-9*y[1];
}

inline double func8LW(double x[], double y[]){  //f(x,y)
  return -9*y[0]+9*y[1];
}

inline int func8CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func8CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
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

inline int func8CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func8CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-6*x[0]+x[1]+x[2]-3*x[3]-9*y[0]-7*y[1]+15;
  constraintValuesListReturn[1]=4*x[1]+5*x[2]+10*x[3]-26;
  constraintValuesListReturn[2]=-9*x[0]+9*x[1]-9*x[2]+5*x[3]-5*y[0]-4*y[1]+5;
  constraintValuesListReturn[3]=5*x[0]+3*x[1]+x[2]+9*x[3]+y[0]+5*y[1]-32;
  constraintValuesListReturn[4]=-y[0];
  constraintValuesListReturn[5]=-y[1];
  return 1;
}


inline int func8CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=-9+(dualNeq[0]*(-9) + dualNeq[1]*(0) + dualNeq[2]*(-5) + dualNeq[3]*(1) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=9+(dualNeq[0]*(-7) + dualNeq[1]*(0) + dualNeq[2]*(-4) + dualNeq[3]*(5) + dualNeq[4]*(0) + dualNeq[5]*(-1));
    return 1;						
}

inline int func8SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound8[12]={0,10e5,0,10e5,0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

/* Função A9 - solucao (0,2.44,10,0,10,8.74,5.25,10,0,10,3.73,10,10,10,0,0)*/

inline double func9UP(double x[], double y[]){  //F(x,y)
  return 12*x[0]-x[1]-12*x[2]+13*x[3]+2*x[5]-5*x[7]+6*x[8]-11*x[9]-5*y[0]-6*y[1]-4*y[2]-7*y[3];
}

inline double func9LW(double x[], double y[]){  //f(x,y)
  return 3*y[0]-2*y[1]-3*y[2]-3*y[3]+y[4]+6*y[5];
}

inline int func9CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int func9CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
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

inline int func9CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int func9CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
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
	constraintValuesListReturn[2*i+2]=-y[i];
	constraintValuesListReturn[2*i+3]=y[i]-10;
  }

  return 1;
}


inline int func9CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
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

inline int func9SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

const double bound9[32]={0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10};  //Bounds x, y

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

const double boundB1[8]={0,50,0,50,-10,20,-10,20};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B2 - solucao (0,30,-10,10)*/

inline double funcB2UP(double x[], double y[]){  //F(x,y)
  return -x[0]*x[0] - 3*x[1]*x[1] - 4*y[0] + y[1]*y[1];
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

const double boundB2[8]={0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y


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

    constraintValuesListReturn[0]=2*y[0]+(dualNeq[0]*(0.4) + dualNeq[1]*(0.6) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(1) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[1]=2*y[1]+(dualNeq[0]*(0.7) + dualNeq[1]*(0.3) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(1) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[2]=2*y[2]+(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(0.4) + dualNeq[3]*(0.6) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(-1) + dualNeq[9]*(1) + dualNeq[10]*(0) + dualNeq[11]*(0));

    constraintValuesListReturn[3]=2*y[3]+(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(0.7) + dualNeq[3]*(0.3) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(0) + dualNeq[9]*(0) + dualNeq[10]*(-1) + dualNeq[11]*(1));
  
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

const double boundB4[18]={0,10e5,0,10e5,0,10e5,0,10e5,0,10e5,0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função B5 - solucao (1.95, 8.05, 0.00, 0.97, 0.97, 1.31, 6.74, 0, 0)*/

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

const double boundB5[8]={0,10e5,0,10e5,0,10e5,0,10e5};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


const Function listFunction[DEFINEfunctionListSize]={{1,1,0,1,0,4,bound1,func1UP,func1LW,func1CTREQUP,func1CTRNEQUP, func1CTREQLW,func1CTRNEQLW,func1CTKKT,func1SimplexTableauKKT,"funcA1"},
						     {1,1,0,0,0,7,bound2,func2UP,func2LW,func2CTREQUP,func2CTRNEQUP, func2CTREQLW,func2CTRNEQLW,func2CTKKT,func2SimplexTableauKKT,"funcA2"},
						     {1,1,0,0,0,8,bound3,func3UP,func3LW,func3CTREQUP,func3CTRNEQUP, func3CTREQLW,func3CTRNEQLW,func3CTKKT,func3SimplexTableauKKT,"funcA3"},
						     {2,3,0,2,0,6,bound4,func4UP,func4LW,func4CTREQUP,func4CTRNEQUP, func4CTREQLW,func4CTRNEQLW,func4CTKKT,func4SimplexTableauKKT,"funcA4"},
						     {1,2,0,2,0,5,bound5,func5UP,func5LW,func5CTREQUP,func5CTRNEQUP, func5CTREQLW,func5CTRNEQLW,func5CTKKT,func5SimplexTableauKKT,"funcA5"},
						     {2,2,0,3,0,4,bound6,func6UP,func6LW,func6CTREQUP,func6CTRNEQUP, func6CTREQLW,func6CTRNEQLW,func6CTKKT,func6SimplexTableauKKT,"funcA6"},
						     {2,3,0,3,0,6,bound7,func7UP,func7LW,func7CTREQUP,func7CTRNEQUP, func7CTREQLW,func7CTRNEQLW,func7CTKKT,func7SimplexTableauKKT,"funcA7"},
						     {4,2,0,10,0,6,bound8,func8UP,func8LW,func8CTREQUP,func8CTRNEQUP, func8CTREQLW,func8CTRNEQLW,func8CTKKT,func8SimplexTableauKKT,"funcA8"},
						     {10,6,0,22,0,19,bound9,func9UP,func9LW,func9CTREQUP,func9CTRNEQUP, func9CTREQLW,func9CTRNEQLW,func9CTKKT,func9SimplexTableauKKT,"funcA9"},
						     {2,2,0,5,0,6,boundB1,funcB1UP,funcB1LW,funcB1CTREQUP,funcB1CTRNEQUP, funcB1CTREQLW,funcB1CTRNEQLW,funcB1CTKKT,funcB1SimplexTableauKKT,"funcB1"},
						     {2,2,0,3,0,4,boundB2,funcB2UP,funcB2LW,funcB2CTREQUP,funcB2CTRNEQUP, funcB2CTREQLW,funcB2CTRNEQLW,funcB2CTKKT,funcB2SimplexTableauKKT,"funcB2"},
						     {4,4,0,9,0,12,boundB3,funcB3UP,funcB3LW,funcB3CTREQUP,funcB3CTRNEQUP, funcB3CTREQLW,funcB3CTRNEQLW,funcB3CTKKT,funcB3SimplexTableauKKT,"funcB3"},
						     {3,6,0,4,0,9,boundB4,funcB4UP,funcB4LW,funcB4CTREQUP,funcB4CTRNEQUP, funcB4CTREQLW,funcB4CTRNEQLW,funcB4CTKKT,funcB4SimplexTableauKKT,"funcB4"},
						     {2,2,0,3,0,4,boundB5,funcB5UP,funcB5LW,funcB5CTREQUP,funcB5CTRNEQUP, funcB5CTREQLW,funcB5CTRNEQLW,funcB5CTKKT,funcB5SimplexTableauKKT,"funcB5"}};



#endif
