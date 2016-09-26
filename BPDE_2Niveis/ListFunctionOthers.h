#ifndef LISTFUNCTIONOTHERS_INCLUDED
#define LISTFUNCTIONOTHERS_INCLUDED  
#define DEFINEfunctionListSizeOthers 21
#include <float.h>
#include <cmath>
#include <stdlib.h>


typedef struct functionPrototypeOthers{
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
} FunctionOthers;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/

//An Evolutionary Algorithm for Solving Nonlinear Bilevel Programming Based on a New Constraint-Handling Scheme

/* Função 5 */

inline double funcO1UP(double x[], double y[]){  //F(x,y)
  double r=0.1;
  
  return r*(x[0]*x[0]+x[1]*x[1])-3*y[0]-4*y[1]+0.5*(y[0]*y[0]+y[1]*y[1]);
}

inline double funcO1LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,-2},{-2,5}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*x[0]+y[1]*x[1]); 
}

inline int funcO1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,-2},{-2,5}};
  
    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,-2},{-2,5}};
  
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1]));
	return 1;
}

inline int funcO1LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,-2},{-2,5}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-x[0];  
      matrixCB[1]=-x[1];
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}


//const double boundO1[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO1[8]={0,1e2,0,1e2,0,3,0,3};  //Bounds x, y




/* Função 6 */

inline double funcO2UP(double x[], double y[]){  //F(x,y)
  double r=1;
  
  return r*(x[0]*x[0]+x[1]*x[1])-3*y[0]-4*y[1]+0.5*(y[0]*y[0]+y[1]*y[1]);
}

inline double funcO2LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,-2},{-2,5}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*x[0]+y[1]*x[1]); 
}

inline int funcO2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,-2},{-2,5}};
  
    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,-2},{-2,5}};
  
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1]));
	return 1;
}

inline int funcO2LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,-2},{-2,5}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-x[0];  
      matrixCB[1]=-x[1];
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}



//const double boundO2[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO2[8]={0,1e2,0,1e2,0,3,0,3};  //Bounds x, y



/* Função 7 */

inline double funcO3UP(double x[], double y[]){  //F(x,y)
  double r=0;
  
  return r*(x[0]*x[0]+x[1]*x[1])-3*y[0]-4*y[1]+0.5*(y[0]*y[0]+y[1]*y[1]);
}

inline double funcO3LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,3},{3,10}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*x[0]+y[1]*x[1]); 
}

inline int funcO3CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO3CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO3CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO3CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO3CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,3},{3,10}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO3SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,3},{3,10}};  
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1]));
	return 1;
}



inline int funcO3LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,3},{3,10}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-x[0];  
      matrixCB[1]=-x[1];
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}



//const double boundO3[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO3[8]={0,1e2,0,1e2,0,3,0,3};  //Bounds x, y




/* Função 8 */

inline double funcO4UP(double x[], double y[]){  //F(x,y)
  double r=0.1;
  
  return r*(x[0]*x[0]+x[1]*x[1])-3*y[0]-4*y[1]+0.5*(y[0]*y[0]+y[1]*y[1]);
}

inline double funcO4LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,3},{3,10}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*x[0]+y[1]*x[1]); 
}

inline int funcO4CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO4CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO4CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO4CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO4CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,3},{3,10}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO4SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,3},{3,10}};  
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[1]));
	return 1;
}


inline int funcO4LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,3},{3,10}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-x[0];  
      matrixCB[1]=-x[1];
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

//const double boundO4[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO4[8]={0,1e2,0,1e2,0,3,0,3};  //Bounds x, y



/* Função 9 */

inline double funcO5UP(double x[], double y[]){  //F(x,y)
  double r=0.1;
  
  return r*(x[0]*x[0]+x[1]*x[1])-3*y[0]-4*y[1]+0.5*(y[0]*y[0]+y[1]*y[1]);
}

inline double funcO5LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,3},{3,10}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(-x[0]+2*x[1])+y[1]*(3*x[0]-3*x[1])); 
}

inline int funcO5CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO5CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO5CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO5CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO5CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,3},{3,10}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (-x[0]+2*x[1])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (3*x[0]-3*x[1])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO5SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,3},{3,10}};  
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (-x[0]+2*x[1]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (3*x[0]-3*x[1]));
	return 1;
}



inline int funcO5LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,3},{3,10}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(-x[0]+2*x[1]);  
      matrixCB[1]=-(3*x[0]-3*x[1]);
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}


//const double boundO5[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO5[8]={0,1e2,0,1e2,0,3,0,3};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

//On the Numerical Solution of a Class of Stackelberg Problems

/*Função 4.1*/

inline double funcO6UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO6LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,0},{0,1}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO6CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO6CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO6CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO6CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO6CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,0},{0,1}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO6SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,0},{0,1}};
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO6LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,0},{0,1}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

//const double boundO6[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO6[6]={0,1e2,0,3,0,3};  //Bounds x, y



/*Função 4.2*/

inline double funcO7UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO7LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1+x[0],0},{0,0}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO7CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO7CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO7CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO7CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO7CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1+x[0],0},{0,0}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO7SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1+x[0],0},{0,0}};
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO7LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1+x[0],0},{0,0}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

//const double boundO7[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO7[6]={0,1e2,0,3,0,3};  //Bounds x, y



/*Função 4.3*/

inline double funcO8UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO8LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1+x[0],0},{0,1+0.1*x[0]}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO8CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO8CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO8CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO8CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-0.333*y[0]+y[1]-2;
  constraintValuesListReturn[1]=y[0]-0.333*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO8CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1+x[0],0},{0,1+0.1*x[0]}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO8SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1+x[0],0},{0,1+0.1*x[0]}};
	
	tableau[0]=-0.333;
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333;
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO8LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1+x[0],0},{0,1+0.1*x[0]}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333;
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333;
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=2;
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

//const double boundO8[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO8[6]={0,1e2,0,3,0,3};  //Bounds x, y



/*Função 4.4*/

inline double funcO9UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO9LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1,0},{0,1}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO9CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO9CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO9CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO9CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=(-0.333+0.1*x[0])*y[0]+y[1]-x[0];
  constraintValuesListReturn[1]=y[0]+(-0.333-0.1*x[0])*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO9CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1,0},{0,1}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333+0.1*x[0]) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333-0.1*x[0]) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO9SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1,0},{0,1}};
	
	tableau[0]=-0.333+0.1*x[0];
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333-0.1*x[0];
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO9LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1,0},{0,1}};					//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333+0.1*x[0];
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333-0.1*x[0];
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=x[0];
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

const double boundO9[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y Sem limitação aparentemente



/*Função 4.5*/

inline double funcO10UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO10LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1+x[0],0},{0,1}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO10CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO10CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO10CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO10CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=(-0.333+0.1*x[0])*y[0]+y[1]-x[0];
  constraintValuesListReturn[1]=y[0]+(-0.333-0.1*x[0])*y[1]-2;
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO10CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1+x[0],0},{0,1}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333+0.1*x[0]) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333-0.1*x[0]) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO10SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1+x[0],0},{0,1}};
	
	tableau[0]=-0.333+0.1*x[0];
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333-0.1*x[0];
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO10LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1+x[0],0},{0,1}};				//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333+0.1*x[0];
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333-0.1*x[0];
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=x[0];
      matrixCB[3]=2;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

const double boundO10[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y  Sem limitação aparentemente



/*Função 4.6*/

inline double funcO11UP(double x[], double y[]){  //F(x,y)  
  return 0.5*((y[0]-3)*(y[0]-3) + (y[1]-4)*(y[1]-4));
}

inline double funcO11LW(double x[], double y[]){  //f(x,y)
  double H[2][2]={{1+0.2*x[0],0},{0,1+0.1*x[0]}};
    
  return 0.5*( (y[0]*H[0][0]+y[1]*H[1][0])*y[0] + (y[0]*H[0][1]+y[1]*H[1][1])*y[1]) - (y[0]*(3+1.333*x[0])+y[1]*(x[0])); 
}

inline int funcO11CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO11CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO11CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO11CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=(-0.333+0.1*x[0])*y[0]+y[1]-2+0.1*x[0];
  constraintValuesListReturn[1]=y[0]+(-0.333-0.1*x[0])*y[1]-2+0.1*x[0];
  constraintValuesListReturn[2]=-y[0];
  constraintValuesListReturn[3]=-y[1];
  return 1;
}


inline int funcO11CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    double H[2][2]={{1+0.2*x[0],0},{0,1+0.1*x[0]}};

    constraintValuesListReturn[0]=(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0])) + (dualNeq[0]*(-0.333+0.1*x[0]) + dualNeq[1]*(1) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(-0.333-0.1*x[0]) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    return 1;						
}

inline int funcO11SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	double H[2][2]={{1+0.2*x[0],0},{0,1+0.1*x[0]}};
	
	tableau[0]=-0.333+0.1*x[0];
	tableau[1]=1;
	tableau[2]=-1;
	tableau[3]=0;
	tableau[4]=-(0.5*(2*y[0]*H[0][0] + y[1]*H[1][0] + y[1]*H[0][1]) - (3+1.333*x[0]));
	
	tableau[5]=1;
	tableau[6]=-0.333-0.1*x[0];
	tableau[7]=0;
	tableau[8]=-1;
	tableau[9]=-(0.5*(y[0]*H[1][0] + y[0]*H[0][1] + 2*y[1]*H[1][1]) - (x[0]));
	return 1;
}



inline int funcO11LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      double H[2][2]={{1+0.2*x[0],0},{0,1+0.1*x[0]}};				//(f)  - (-A  0 )(l1)= (-b)
     
      matrixQ[0]=H[0][0];
      matrixQ[1]=H[0][1];
      matrixQ[2]=H[1][0];
      matrixQ[3]=H[1][1];
      
      matrixA[0]=-0.333+0.1*x[0];
      matrixA[1]=1;
      
      matrixA[2]=1;
      matrixA[3]=-0.333-0.1*x[0];
      
      matrixA[4]=-1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      
      matrixCB[0]=-(3+1.333*x[0]);  
      matrixCB[1]=-(x[0]);
      matrixCB[2]=2-0.1*x[0];
      matrixCB[3]=2-0.1*x[0];
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      
      return 1;
}

//const double boundO11[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO11[6]={0,20,0,4,0,4};  //Bounds x, y Aproximadamente


/* --------------------------------------------------------------------------------------------------------------------------------------*/



//Solving convex quadratic bilevel programming problems using an enumeration

/*Função 2.1 Semelhante ao problema J9*/

inline double funcO12UP(double x[], double y[]){  //F(x,y)  
  return (x[0]-1)*(x[0]-1)+2*y[0]*y[0]-2*x[0];
}

inline double funcO12LW(double x[], double y[]){  //f(x,y)    
  return (2*y[0]-x[0])*(2*y[0]-x[0])+(2*y[1]-1)*(2*y[1]-1)+x[0]*y[0];
}

inline int funcO12CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO12CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO12CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO12CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=4*x[0]+5*y[0]+4*y[1]-12;
  constraintValuesListReturn[1]=-4*x[0]-5*y[0]+4*y[1]+4;
  constraintValuesListReturn[2]=4*x[0]-4*y[0]+5*y[1]-4;
  constraintValuesListReturn[3]=-4*x[0]+4*y[0]+5*y[1]-4;
  constraintValuesListReturn[4]=-x[0];
  constraintValuesListReturn[5]=-y[0];
  constraintValuesListReturn[6]=-y[1];
  
  
  return 1;
}


inline int funcO12CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
  
    constraintValuesListReturn[0]=(2*2*(2*y[0]-x[0])+x[0]) + (dualNeq[0]*(5) + dualNeq[1]*(-5) + dualNeq[2]*(-4) + dualNeq[3]*(4) + dualNeq[4]*(0) + dualNeq[5]*(-1) + dualNeq[6]*(0));
    constraintValuesListReturn[1]=(2*2*(2*y[1]-1)) + (dualNeq[0]*(4) + dualNeq[1]*(4) + dualNeq[2]*(5) + dualNeq[3]*(5) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1));
    return 1;						
}

inline int funcO12SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))	
	tableau[0]=5;
	tableau[1]=-5;
	tableau[2]=-4;
	tableau[3]=4;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=0;
	tableau[7]=-(2*2*(2*y[0]-x[0])+x[0]);
	
	tableau[8]=4;
	tableau[9]=4;
	tableau[10]=5;
	tableau[11]=5;
	tableau[12]=0;
	tableau[13]=0;
	tableau[14]=-1;
	tableau[15]=-(2*2*(2*y[1]-1));
	
	return 1;
}



inline int funcO12LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
								//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[0]=8;
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=8;
      
      matrixA[0]=5;
      matrixA[1]=4;
      
      matrixA[2]=-5;
      matrixA[3]=4;
      
      matrixA[4]=-4;
      matrixA[5]=5;
      
      matrixA[6]=4;
      matrixA[7]=5;

      matrixA[8]=0;
      matrixA[9]=0;
      
      matrixA[10]=-1;
      matrixA[11]=0;
      
      matrixA[12]=0;
      matrixA[13]=-1;
      
      matrixCB[0]=-4*x[0]+x[0];  
      matrixCB[1]=-4;
      matrixCB[2]=-4*x[0]+12;
      matrixCB[3]=4*x[0]-4;
      matrixCB[4]=-4*x[0]+4;
      matrixCB[5]=4*x[0]+4;
      matrixCB[6]=0;
      matrixCB[7]=0;
      matrixCB[8]=0;     
      
      return 1;
}

//const double boundO12[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundO12[6]={0,2,0,2,0,1};  //Bounds x, y



/*Função 5.1*/

inline double funcO13UP(double x[], double y[]){  //F(x,y)  
  return (x[0]-30)*(x[0]-30)+(x[1]-20)*(x[1]-20)-20*y[0]+20*y[1];
}

inline double funcO13LW(double x[], double y[]){  //f(x,y)    
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]);
}

inline int funcO13CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO13CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+2*x[1]-30;
  constraintValuesListReturn[1]=-(x[0]+x[1]-20);
  constraintValuesListReturn[2]=-x[0];
  constraintValuesListReturn[3]=x[0]-15;
  constraintValuesListReturn[4]=-x[1];
  constraintValuesListReturn[5]=x[1]-15;
  
  return 1;
}

inline int funcO13CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO13CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0];
  constraintValuesListReturn[1]=y[0]-15;
  constraintValuesListReturn[2]=-y[1];
  constraintValuesListReturn[3]=y[1]-15;
  
  return 1;
}


inline int funcO13CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
  
    constraintValuesListReturn[0]=(-1*2*(x[0]-y[0])) + (dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(-1*2*(x[1]-y[1])) + (dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1));
    return 1;						
}

inline int funcO13SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))	
	tableau[0]=-1;
	tableau[1]=1;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=-(-1*2*(x[0]-y[0]));
	
	tableau[5]=0;
	tableau[6]=0;
	tableau[7]=-1;
	tableau[8]=1;
	tableau[9]=-(-1*2*(x[1]-y[1]));
	
	return 1;
}



inline int funcO13LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
								//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[0]=2;
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
      matrixCB[3]=15;
      matrixCB[4]=0;
      matrixCB[5]=15;    
      
      return 1;
}

const double boundO13[8]={0,15,0,15,0,15,0,15};  //Bounds x, y



/*Função 5.2 PROBLEMA NA FO UP*/

inline double funcO14UP(double x[], double y[]){  //F(x,y)  
  return y[0]*y[0]+y[1]*y[1]-y[0]*y[2]-4*y[1]-7*x[0]+4*x[1];
}

inline double funcO14LW(double x[], double y[]){  //f(x,y)    
  return y[0]*y[0]+0.5*y[1]*y[1]+0.5*y[2]*y[2]+y[0]*y[1]+(1-3*x[0])*y[0] + (1+x[1])*y[1];
}

inline int funcO14CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO14CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]-1;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-x[1];
  
  return 1;
}

inline int funcO14CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO14CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=2*y[0]+y[1]-y[2]+x[0]-2*x[1]+2;
  constraintValuesListReturn[1]=-y[0];
  constraintValuesListReturn[2]=-y[1];
  constraintValuesListReturn[3]=-y[2];
  
  return 1;
}


inline int funcO14CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
  
    constraintValuesListReturn[0]=(2*y[0]+y[1]+(1-3*x[0])) + (dualNeq[0]*(2) + dualNeq[1]*(-1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(y[1]+y[0]+(1+x[0])) + (dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(0));
    constraintValuesListReturn[2]=(y[2]) + (dualNeq[0]*(-1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(-1));
    
    return 1;						
}

inline int funcO14SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))	
	tableau[0]=2;
	tableau[1]=-1;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=-(2*y[0]+y[1]+(1-3*x[0]));
	
	tableau[5]=1;
	tableau[6]=0;
	tableau[7]=-1;
	tableau[8]=0;
	tableau[9]=-(y[1]+y[0]+(1+x[0]));
	
	tableau[10]=-1;
	tableau[11]=0;
	tableau[12]=0;
	tableau[13]=-1;
	tableau[14]=-(y[2]);
	
	
	return 1;
}



inline int funcO14LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
								//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[0]=2;
      matrixQ[1]=1;
      matrixQ[2]=0;
      matrixQ[3]=1;
      matrixQ[4]=1;
      matrixQ[5]=0;
      matrixQ[6]=0;
      matrixQ[7]=0;
      matrixQ[8]=1;
      
      
      matrixA[0]=2;
      matrixA[1]=1;
      matrixA[2]=-1;
      
      matrixA[3]=-1;
      matrixA[4]=0;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;
      matrixA[8]=0;
      
     
      matrixA[9]=0;
      matrixA[10]=0;
      matrixA[11]=-1;

      matrixCB[0]=(1-3*x[0]);  
      matrixCB[1]=(1+x[0]);
      matrixCB[2]=0;
      matrixCB[3]=-x[0]+2*x[1]-2;
      matrixCB[4]=0;
      matrixCB[5]=0;    
      matrixCB[6]=0;    
      
      return 1;
}

const double boundO14[10]={0,1,0,1,0,1e2,0,1e2,0,1e2};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/

//BIPA Arquivo Instâncias 

/*DEMPE92   PROBLEMA NA RESTRIÇÂO*/

inline double funcO15UP(double x[], double y[]){  //F(x,y)  
  return (x[0]-3.5)*(x[0]-3.5) + (y[0]-4)*(y[0]-4);
}

inline double funcO15LW(double x[], double y[]){  //f(x,y)    
  return (y[0]-3)*(y[0]-3);
}

inline int funcO15CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO15CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcO15CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO15CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=y[0]*y[0]-x[0];
  
  return 1;
}


inline int funcO15CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(y[0]-3)) + (dualNeq[0]*(2*y[0]));
    
    return 1;						
}

inline int funcO15SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))	
	tableau[0]=2*y[0];
	tableau[1]=-(2*(y[0]-3));
	
	return 1;
}



inline int funcO15LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[0]=2;  
      
      matrixA[0]=2;//*y[0];

      matrixCB[0]=-6;  
      matrixCB[1]=x[0];
      
      return 1;
}

const double boundO15[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y  FUNÇÂO COM PROBLEMAS


/* --------------------------------------------------------------------------------------------------------------------------------------*/

//A Hybrid Genetic Algorithm for Solving a Class of Nonlinear Bilevel Programming Problems

/* Função F19*/

inline double funcO16UP(double x[], double y[]){  //F(x,y)
  return fabs(sin((x[0]-30)*(x[0]-30)+(x[1]-20)*(x[1]-20)-20*y[0]+20*y[1]-225));
}

inline double funcO16LW(double x[], double y[]){  //f(x,y)
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]);
}

inline int funcO16CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO16CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=- (x[0]+2*x[1]-30);
  constraintValuesListReturn[1]=x[0]+x[1]-25;
  constraintValuesListReturn[2]=x[1]-15;

  return 1;
}

inline int funcO16CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO16CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0];
  constraintValuesListReturn[1]=y[0]-10;
  constraintValuesListReturn[2]=-y[1];
  constraintValuesListReturn[3]=y[1]-10;
  return 1;
}


inline int funcO16CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=(-1*2*(x[0]-y[0])) +(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(-1*2*(x[1]-y[1])) +(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1));
    return 1;						
}

inline int funcO16SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=-1;
	tableau[1]=1;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=-(-1*2*(x[0]-y[0]));
	
	tableau[5]=0;
	tableau[6]=0;
	tableau[7]=-1;
	tableau[8]=1;
	tableau[9]=-(-1*2*(x[1]-y[1]));
	return 1;
}

inline int funcO16LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)
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

//const double boundO16[8]={-1e2,1e2,-1e2,15,0,10,0,10};  //Bounds x, y
const double boundO16[8]={0,20,5,15,0,10,0,10};  //Bounds x, y




/* Função F20*/

inline double funcO17UP(double x[], double y[]){  //F(x,y) 
  return fabs(sin(2*x[0]+2*x[1]-3*y[0]-3*y[1]-60));
}

inline double funcO17LW(double x[], double y[]){  //f(x,y) 
  return (y[0]-x[0]+20)*(y[0]-x[0]+20) + (y[1]-x[1]+20)*(y[1]-x[1]+20);
}

inline int funcO17CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO17CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]+y[0]-2*y[1]-40;

  return 1;
}

inline int funcO17CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO17CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=2*y[0]-x[0]+10;
  constraintValuesListReturn[1]=2*y[1]-x[1]+10;
  constraintValuesListReturn[2]=-x[0];
  constraintValuesListReturn[3]=x[0]-50;
  constraintValuesListReturn[4]=-x[1];
  constraintValuesListReturn[5]=x[1]-50;
  constraintValuesListReturn[6]=-y[0]-10;
  constraintValuesListReturn[7]=y[0]-20;
  constraintValuesListReturn[8]=-y[1]-10;
  constraintValuesListReturn[9]=y[1]-20;
  
  return 1;
}


inline int funcO17CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) 
    constraintValuesListReturn[0]=(2*(y[0]-x[0]+20)) + (dualNeq[0]*(2) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(1) + dualNeq[8]*(0) + dualNeq[9]*(0));
    constraintValuesListReturn[1]=(2*(y[1]-x[1]+20)) + (dualNeq[0]*(0) + dualNeq[1]*(2) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(-1) + dualNeq[9]*(1));

    return 1;						
}

inline int funcO17SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=2;
	tableau[1]=0;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=0;
	tableau[6]=-1;
	tableau[7]=1;
	tableau[8]=0;
	tableau[9]=0;
	tableau[10]=-(2*(y[0]-x[0]+20)) ;
	
	tableau[11]=0;
	tableau[12]=2;
	tableau[13]=0;
	tableau[14]=0;
	tableau[15]=0;
	tableau[16]=0;
	tableau[17]=0;
	tableau[18]=0;
	tableau[19]=-1;
	tableau[20]=1;
	tableau[21]=-(2*(y[1]-x[1]+20));
	
	return 1;
}


/*LEMKE COM PROBLMA POR y NEGATIVO*/
inline int funcO17LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=2;
      
      matrixA[0]=2;
      matrixA[1]=0;
      
      matrixA[2]=0;
      matrixA[3]=2;
      
      matrixA[4]=0;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=0;

      matrixA[8]=0;
      matrixA[9]=0;

      matrixA[10]=0;
      matrixA[11]=0;

      matrixA[12]=-1;
      matrixA[13]=0;

      matrixA[14]=1;
      matrixA[15]=0;

      matrixA[16]=0;
      matrixA[17]=-1;         
      
      matrixA[16]=0;
      matrixA[17]=1; 
      
      matrixCB[0]=-2*x[0]+40;  
      matrixCB[1]=-2*x[1]+40;
      matrixCB[2]=x[0]-10;
      matrixCB[3]=x[1]-10;
      matrixCB[4]=0;
      matrixCB[5]=50;
      matrixCB[6]=0;
      matrixCB[7]=50;
      matrixCB[8]=10;
      matrixCB[9]=20;
      matrixCB[10]=10;
      matrixCB[11]=20;
      
      return 1;
}


const double boundO17[8]={0,50,0,50,-10,20,-10,20};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

//A dual temperature simulated annealing approach for solving bilevel programming problems

/* Função Example 1*/


inline double funcO18UP(double x[], double y[]){  //F(x,y) MAX
  return -(x[0]+3*y[0]-2*y[1]);
}

inline double funcO18LW(double x[], double y[]){  //f(x,y) MAX
  return -(y[0]);
}

inline int funcO18CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO18CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-8;

  return 1;
}

inline int funcO18CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO18CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]-4*y[1]+16);
  constraintValuesListReturn[1]=-(-8*x[0]-3*y[0]+2*y[1]+48);
  constraintValuesListReturn[2]=-(2*x[0]-y[0]+3*y[1]-12);
  constraintValuesListReturn[3]=-y[0];
  constraintValuesListReturn[4]=y[0]-4;

  
  return 1;
}


inline int funcO18CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=(1) - (dualNeq[0]*(1) + dualNeq[1]*(3) + dualNeq[2]*(1) + dualNeq[3]*(-1) + dualNeq[4]*(1));
    constraintValuesListReturn[1]=(0) - (dualNeq[0]*(4) + dualNeq[1]*(-2) + dualNeq[2]*(-3) + dualNeq[3]*(0) + dualNeq[4]*(0));

    return 1;						
}

inline int funcO18SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=1;
	tableau[1]=3;
	tableau[2]=1;
	tableau[3]=-1;
	tableau[4]=1;
	tableau[5]=1;
	
	tableau[6]=4;
	tableau[7]=-2;
	tableau[8]=-3;
	tableau[9]=0;
	tableau[10]=0;
	tableau[11]=0;

	
	return 1;
}



inline int funcO18LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      
      matrixA[0]=1;
      matrixA[1]=4;
      
      matrixA[2]=3;
      matrixA[3]=-2;
      
      matrixA[4]=1;
      matrixA[5]=-3;
      
      matrixA[6]=-1;
      matrixA[7]=0;

      matrixA[8]=1;
      matrixA[9]=0;

   
      matrixCB[0]=-1;  
      matrixCB[1]=0;
      matrixCB[2]=-(-2*x[0]-16);
      matrixCB[3]=-(8*x[0]-48);
      matrixCB[4]=-(-2*x[0]+12);
      matrixCB[5]=0;
      matrixCB[6]=4;
        
      return 1;
}


//const double boundO18[8]={0,8,0,4,-1e2,1e2};  //Bounds x, y
const double boundO18[8]={0,8,0,4,0,8};  //Bounds x, y


/* Função Example 2*/


inline double funcO19UP(double x[], double y[]){  //F(x,y) MAX
  return -(-((x[0]-3)*(x[0]-3) + (y[0]-2)*(y[0]-2)));
}

inline double funcO19LW(double x[], double y[]){  //f(x,y) MAX
  return -(-((y[0]-5)*(y[0]-5)));
}

inline int funcO19CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO19CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-8;

  return 1;
}

inline int funcO19CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO19CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+1);
  constraintValuesListReturn[1]=-(-x[0]+2*y[0]-2);
  constraintValuesListReturn[2]=-(-x[0]-2*y[0]+14);
  
  return 1;
}


inline int funcO19CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=(-2*(y[0]-5)) - (dualNeq[0]*(1) + dualNeq[1]*(-2) + dualNeq[2]*(2));

    return 1;						
}

inline int funcO19SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=1;
	tableau[1]=-2;
	tableau[2]=2;
	tableau[3]=-2*(y[0]-5);
		
	return 1;
}



inline int funcO19LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)

      
      matrixA[0]=1;
      
      matrixA[1]=-2;
      
      matrixA[2]=2;
      
  
      matrixCB[0]=-10;  
      matrixCB[1]=-(-2*x[0]-1);
      matrixCB[2]=-(x[0]+2);
      matrixCB[3]=-(x[0]-14);
        
      return 1;
}


//const double boundO19[8]={0,8,-1e2,1e2};  //Bounds x, y
const double boundO19[8]={0,8,1,5.8};  //Bounds x, y




/* Função Example 3*/


inline double funcO20UP(double x[], double y[]){  //F(x,y) MAX
  return -(-((x[0]-3)*(x[0]-3) + (y[0]-2)*(y[0]-2)));
}

inline double funcO20LW(double x[], double y[]){  //f(x,y) MAX
  return -(-((y[0]-5)*(y[0]-5)));
}

inline int funcO20CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO20CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+1);
  constraintValuesListReturn[1]=-(-x[0]+2*y[0]-2);
  constraintValuesListReturn[2]=-(-x[0]-2*y[0]+14);
  constraintValuesListReturn[3]=-x[0];
  constraintValuesListReturn[4]=x[0]-8;

  return 1;
}

inline int funcO20CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO20CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}


inline int funcO20CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=(-2*(y[0]-5));

    return 1;						
}

inline int funcO20SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=-2*(y[0]-5);
	
	return 1;
}



inline int funcO20LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)

      matrixCB[0]=-10;  
        
      return 1;
}


//const double boundO20[8]={0,8,-1e2,1e2};  //Bounds x, y
const double boundO20[8]={0,8,1,5.8};  //Bounds x, y



/* Função Example 4*/


inline double funcO21UP(double x[], double y[]){  //F(x,y) MAX
  return -(  (-2/5*x[0]*x[0]*x[1]+4*x[1]*x[1])*(y[0]*y[1]) + (-x[1]*x[1]*x[1]+3*x[0]*x[0]*x[1])*(1-y[0])*y[1] + (2*x[1]*x[1]-x[0])*(1-y[1])  );
}

inline double funcO21LW(double x[], double y[]){  //f(x,y) MAX
  return -( (x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1])*(y[0]*y[1]) + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*(1-y[0])*y[1] + 8*x[0]*y[0]*(1-y[1]) );
}

inline int funcO21CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcO21CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-10;
  constraintValuesListReturn[2]=-x[1];
  constraintValuesListReturn[3]=x[1]-10;

  return 1;
}

inline int funcO21CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcO21CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(y[0]+y[1]-1);
  constraintValuesListReturn[1]=-y[0];
  constraintValuesListReturn[2]=y[0]-1;
  constraintValuesListReturn[3]=-y[1];
  constraintValuesListReturn[4]=y[0]-1;

  
  return 1;
}


inline int funcO21CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=( (x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1])*y[1] + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*y[1]*(-1) + 8*x[0]*(1-y[1]) ) - (dualNeq[0]*(-1) + dualNeq[1]*(-1) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(0));
    constraintValuesListReturn[1]=( (x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1])*y[0] + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*(1-y[0]) + 8*x[0]*y[0]*(-1) ) - (dualNeq[0]*(-1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(-1) + dualNeq[4]*(1));

    return 1;						
}

inline int funcO21SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=-1;
	tableau[1]=-1;
	tableau[2]=1;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=(x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1])*y[1] + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*y[1]*(-1) + 8*x[0]*(1-y[1]);
	
	tableau[6]=-1;
	tableau[7]=0;
	tableau[8]=0;
	tableau[9]=-1;
	tableau[10]=1;
	tableau[11]=(x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1])*y[0] + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*(1-y[0]) + 8*x[0]*y[0]*(-1);

	
	return 1;
}



inline int funcO21LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[1]=- ((x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1]) + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*(-1) + 8*x[0]*(-1));
      matrixQ[2]=- ((x[0]*x[0]*x[1]*x[1] + 8*x[1]*x[1]*x[1] - 14*x[0]*x[0] - 5*x[1]) + (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1])*(-1) + 8*x[0]*(-1)) ;
      matrixQ[3]=0;
      
      matrixA[0]=-1;
      matrixA[1]=-1;
      
      matrixA[2]=-1;
      matrixA[3]=0;
      
      matrixA[4]=1;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=-1;

      matrixA[8]=0;
      matrixA[9]=1;

   
      matrixCB[0]=- 8*x[0];  
      matrixCB[1]=- (-x[0]*x[1]*x[1] + 5*x[0]*x[1] + 4*x[1]);
      matrixCB[2]=-1;
      matrixCB[3]=0;
      matrixCB[4]=1;
      matrixCB[5]=0;
      matrixCB[6]=1;
        
      return 1;
}


const double boundO21[8]={0,10,0,10,0,1,0,1};  //Bounds x, y



/* --------------------------------------------------------------------------------------------------------------------------------------*/



const FunctionOthers listFunctionOthers[DEFINEfunctionListSizeOthers]={{2,2,0,0,0,4,boundO1,funcO1UP,funcO1LW,funcO1CTREQUP,funcO1CTRNEQUP, funcO1CTREQLW,funcO1CTRNEQLW,funcO1CTKKT,funcO1SimplexTableauKKT,funcO1LemkeMatrix,"funcO1"},
					           {2,2,0,0,0,4,boundO2,funcO2UP,funcO2LW,funcO2CTREQUP,funcO2CTRNEQUP, funcO2CTREQLW,funcO2CTRNEQLW,funcO2CTKKT,funcO2SimplexTableauKKT,funcO2LemkeMatrix,"funcO2"},
					           {2,2,0,0,0,4,boundO3,funcO3UP,funcO3LW,funcO3CTREQUP,funcO3CTRNEQUP, funcO3CTREQLW,funcO3CTRNEQLW,funcO3CTKKT,funcO3SimplexTableauKKT,funcO3LemkeMatrix,"funcO3"},
					           {2,2,0,0,0,4,boundO4,funcO4UP,funcO4LW,funcO4CTREQUP,funcO4CTRNEQUP, funcO4CTREQLW,funcO4CTRNEQLW,funcO4CTKKT,funcO4SimplexTableauKKT,funcO4LemkeMatrix,"funcO4"},
					           {2,2,0,0,0,4,boundO5,funcO5UP,funcO5LW,funcO5CTREQUP,funcO5CTRNEQUP, funcO5CTREQLW,funcO5CTRNEQLW,funcO5CTKKT,funcO5SimplexTableauKKT,funcO5LemkeMatrix,"funcO5"},
					           {1,2,0,0,0,4,boundO6,funcO6UP,funcO6LW,funcO6CTREQUP,funcO6CTRNEQUP, funcO6CTREQLW,funcO6CTRNEQLW,funcO6CTKKT,funcO6SimplexTableauKKT,funcO6LemkeMatrix,"funcO6"},
					           {1,2,0,0,0,4,boundO7,funcO7UP,funcO7LW,funcO7CTREQUP,funcO7CTRNEQUP, funcO7CTREQLW,funcO7CTRNEQLW,funcO7CTKKT,funcO7SimplexTableauKKT,funcO7LemkeMatrix,"funcO7"},
					           {1,2,0,0,0,4,boundO8,funcO8UP,funcO8LW,funcO8CTREQUP,funcO8CTRNEQUP, funcO8CTREQLW,funcO8CTRNEQLW,funcO8CTKKT,funcO8SimplexTableauKKT,funcO8LemkeMatrix,"funcO8"},
					           {1,2,0,0,0,4,boundO9,funcO9UP,funcO9LW,funcO9CTREQUP,funcO9CTRNEQUP, funcO9CTREQLW,funcO9CTRNEQLW,funcO9CTKKT,funcO9SimplexTableauKKT,funcO9LemkeMatrix,"funcO9"},
					           {1,2,0,0,0,4,boundO10,funcO10UP,funcO10LW,funcO10CTREQUP,funcO10CTRNEQUP, funcO10CTREQLW,funcO10CTRNEQLW,funcO10CTKKT,funcO10SimplexTableauKKT,funcO10LemkeMatrix,"funcO10"},
					           {1,2,0,0,0,4,boundO11,funcO11UP,funcO11LW,funcO11CTREQUP,funcO11CTRNEQUP, funcO11CTREQLW,funcO11CTRNEQLW,funcO11CTKKT,funcO11SimplexTableauKKT,funcO11LemkeMatrix,"funcO11"},
					           {1,2,0,0,0,7,boundO12,funcO12UP,funcO12LW,funcO12CTREQUP,funcO12CTRNEQUP, funcO12CTREQLW,funcO12CTRNEQLW,funcO12CTKKT,funcO12SimplexTableauKKT,funcO12LemkeMatrix,"funcO12"},
					           {2,2,0,6,0,4,boundO13,funcO13UP,funcO13LW,funcO13CTREQUP,funcO13CTRNEQUP, funcO13CTREQLW,funcO13CTRNEQLW,funcO13CTKKT,funcO13SimplexTableauKKT,funcO13LemkeMatrix,"funcO13"},
					           {2,3,0,3,0,4,boundO14,funcO14UP,funcO14LW,funcO14CTREQUP,funcO14CTRNEQUP, funcO14CTREQLW,funcO14CTRNEQLW,funcO14CTKKT,funcO14SimplexTableauKKT,funcO14LemkeMatrix,"funcO14"},
					           {1,1,0,0,0,1,boundO15,funcO15UP,funcO15LW,funcO15CTREQUP,funcO15CTRNEQUP, funcO15CTREQLW,funcO15CTRNEQLW,funcO15CTKKT,funcO15SimplexTableauKKT,funcO15LemkeMatrix,"funcO15"},
					           {2,2,0,3,0,4,boundO16,funcO16UP,funcO16LW,funcO16CTREQUP,funcO16CTRNEQUP, funcO16CTREQLW,funcO16CTRNEQLW,funcO16CTKKT,funcO16SimplexTableauKKT,funcO16LemkeMatrix,"funcO16"},
					           {2,2,0,1,0,10,boundO17,funcO17UP,funcO17LW,funcO17CTREQUP,funcO17CTRNEQUP, funcO17CTREQLW,funcO17CTRNEQLW,funcO17CTKKT,funcO17SimplexTableauKKT,funcO17LemkeMatrix,"funcO17"},
					           {1,2,0,2,0,5,boundO18,funcO18UP,funcO18LW,funcO18CTREQUP,funcO18CTRNEQUP, funcO18CTREQLW,funcO18CTRNEQLW,funcO18CTKKT,funcO18SimplexTableauKKT,funcO18LemkeMatrix,"funcO18"},
					           {1,1,0,2,0,3,boundO19,funcO19UP,funcO19LW,funcO19CTREQUP,funcO19CTRNEQUP, funcO19CTREQLW,funcO19CTRNEQLW,funcO19CTKKT,funcO19SimplexTableauKKT,funcO19LemkeMatrix,"funcO19"},
					           {1,1,0,5,0,0,boundO20,funcO20UP,funcO20LW,funcO20CTREQUP,funcO20CTRNEQUP, funcO20CTREQLW,funcO20CTRNEQLW,funcO20CTKKT,funcO20SimplexTableauKKT,funcO20LemkeMatrix,"funcO20"},
					           {2,2,0,4,0,5,boundO21,funcO21UP,funcO21LW,funcO21CTREQUP,funcO21CTRNEQUP, funcO21CTREQLW,funcO21CTRNEQLW,funcO21CTKKT,funcO21SimplexTableauKKT,funcO21LemkeMatrix,"funcO21"}
						  
};



#endif
