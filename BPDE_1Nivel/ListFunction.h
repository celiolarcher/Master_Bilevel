#ifndef LISTFUNCTION_INCLUDED
#define LISTFUNCTION_INCLUDED  
#define DEFINEfunctionListSize 9
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
    int (*funcSimplexTableauKKT)(double x[], double y[],double tableau[]);
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


/* Função 3 - solucao (17.45,10.91)*/

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


/* Função 4 - solucao (0,0.9,0,0.6,0.4)*/

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


/* Função 5 - solucao (1,9,0)*/

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


/* Função 6 - solucao (2,0,1.5,0)*/

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


/* Função 7 - solucao (0.5,0.8,0,0.2,0.8)*/

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


/* Função 8 - solucao (1.55,0.78,0.16,2.21,1.89,0)*/

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

/* Função 9 - solucao (0,2.44,10,0,10,8.74,5.25,10,0,10,3.73,10,10,10,0,0)*/

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


const Function listFunction[DEFINEfunctionListSize]={{1,1,0,1,0,4,bound1,func1UP,func1LW,func1CTREQUP,func1CTRNEQUP, func1CTREQLW,func1CTRNEQLW,func1CTKKT,func1SimplexTableauKKT,"func1"},
						     {1,1,0,0,0,7,bound2,func2UP,func2LW,func2CTREQUP,func2CTRNEQUP, func2CTREQLW,func2CTRNEQLW,func2CTKKT,func2SimplexTableauKKT,"func2"},
						     {1,1,0,0,0,8,bound3,func3UP,func3LW,func3CTREQUP,func3CTRNEQUP, func3CTREQLW,func3CTRNEQLW,func3CTKKT,func3SimplexTableauKKT,"func3"},
						     {2,3,0,2,0,6,bound4,func4UP,func4LW,func4CTREQUP,func4CTRNEQUP, func4CTREQLW,func4CTRNEQLW,func4CTKKT,func4SimplexTableauKKT,"func4"},
						     {1,2,0,2,0,5,bound5,func5UP,func5LW,func5CTREQUP,func5CTRNEQUP, func5CTREQLW,func5CTRNEQLW,func5CTKKT,func5SimplexTableauKKT,"func5"},
						     {2,2,0,3,0,4,bound6,func6UP,func6LW,func6CTREQUP,func6CTRNEQUP, func6CTREQLW,func6CTRNEQLW,func6CTKKT,func6SimplexTableauKKT,"func6"},
						     {2,3,0,3,0,6,bound7,func7UP,func7LW,func7CTREQUP,func7CTRNEQUP, func7CTREQLW,func7CTRNEQLW,func7CTKKT,func7SimplexTableauKKT,"func7"},
						     {4,2,0,10,0,6,bound8,func8UP,func8LW,func8CTREQUP,func8CTRNEQUP, func8CTREQLW,func8CTRNEQLW,func8CTKKT,func8SimplexTableauKKT,"func8"},
						     {10,6,0,22,0,19,bound9,func9UP,func9LW,func9CTREQUP,func9CTRNEQUP, func9CTREQLW,func9CTRNEQLW,func9CTKKT,func9SimplexTableauKKT,"func9"}};



#endif
