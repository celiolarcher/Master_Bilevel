#ifndef LISTFUNCTIONJAQ_INCLUDED
#define LISTFUNCTIONJAQ_INCLUDED  
#define DEFINEfunctionListSizeJaq 18
#include <float.h>
#include <cmath>
#include <stdlib.h>


typedef struct functionPrototypeJAQ{
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
} FunctionJAQ;

//Function listFunction[DEFINEfunctionListSize];

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J1 - solucao (10,10)*/

inline double funcJ1UP(double x[], double y[]){  //F(x,y)
  return x[0]*x[0]+(y[0]-10)*(y[0]-10);
}

inline double funcJ1LW(double x[], double y[]){  //f(x,y)
  return (x[0]+2*y[0]-30)*(x[0]+2*y[0]-30);
}

inline int funcJ1CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ1CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]+y[0];
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=x[0]-15;

  return 1;
}

inline int funcJ1CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ1CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+y[0]-20;
  constraintValuesListReturn[1]=-y[0];
  constraintValuesListReturn[2]=y[0]-20;
  return 1;
}


inline int funcJ1CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=(2*2*(x[0]+2*y[0]-30)) +(dualNeq[0]*(1) + dualNeq[1]*(-1) + dualNeq[2]*(1));
    return 1;						
}

inline int funcJ1SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
	tableau[0]=1;
	tableau[1]=-1;
	tableau[2]=1;
	tableau[3]=-(2*2*(x[0]+2*y[0]-30));
	return 1;
}

inline int funcJ1LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=8;							//(f)  - (-A  0 )(l1)= (-b)
      
      matrixA[0]=1;
      
      matrixA[1]=-1;
      
      matrixA[2]=1;
           
      matrixCB[0]=4*x[0]-120;  
      matrixCB[1]=20-x[0];
      matrixCB[2]=0;
      matrixCB[3]=20;
      
      return 1;
}

const double boundJ1[4]={0,15,0,20};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J2 - solucao (10,10)*/

inline double funcJ2UP(double x[], double y[]){  //F(x,y)
  return (x[0]-30)*(x[0]-30)+(x[1]-20)*(x[1]-20)-20*y[0]+20*y[1];
}

inline double funcJ2LW(double x[], double y[]){  //f(x,y)
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]);
}

inline int funcJ2CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ2CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=- (x[0]+2*x[1]-30);
  constraintValuesListReturn[1]=x[0]+x[1]-25;
  constraintValuesListReturn[2]=x[1]-15;

  return 1;
}

inline int funcJ2CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ2CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0];
  constraintValuesListReturn[1]=y[0]-10;
  constraintValuesListReturn[2]=-y[1];
  constraintValuesListReturn[3]=y[1]-10;
  return 1;
}


inline int funcJ2CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)
    constraintValuesListReturn[0]=(-1*2*(x[0]-y[0])) +(dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0));
    constraintValuesListReturn[1]=(-1*2*(x[1]-y[1])) +(dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1));
    return 1;						
}

inline int funcJ2SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))
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

inline int funcJ2LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
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

//const double boundJ2[8]={-1e2,1e2,-1e2,15,0,10,0,10};  //Bounds x, y
const double boundJ2[8]={0,20,5,15,0,10,0,10};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

/* Função J3 - solucao (10,10)*/

inline double funcJ3UP(double x[], double y[]){  //F(x,y) MAX
  return - (8*x[0]+4*x[1]-4*y[0]+40*y[1]+4*y[2]);
}

inline double funcJ3LW(double x[], double y[]){  //f(x,y) MAX
  return -(-x[0]-2*x[1]-y[0]-y[1]-2*y[2]);
}

inline int funcJ3CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ3CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcJ3CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ3CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-y[0]+y[1]+y[2]-1;
  constraintValuesListReturn[1]=2*x[0]-y[0]+2*y[1]-0.5*y[2]-1;
  constraintValuesListReturn[2]=2*x[1]+2*y[0]-y[1]-0.5*y[2]-1;
  constraintValuesListReturn[3]=-x[0];
  constraintValuesListReturn[4]=-x[1];
  constraintValuesListReturn[5]=-y[0];
  constraintValuesListReturn[6]=-y[1];
  constraintValuesListReturn[7]=-y[2];
  
  return 1;
}


inline int funcJ3CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=(-1) - (dualNeq[0]*(-1) + dualNeq[1]*(-1) + dualNeq[2]*(2) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1) + dualNeq[6]*(0) + dualNeq[7]*(0));
    constraintValuesListReturn[1]=(-1) - (dualNeq[0]*(1) + dualNeq[1]*(2) + dualNeq[2]*(-1) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(0));
    constraintValuesListReturn[2]=(-2) - (dualNeq[0]*(1) + dualNeq[1]*(-0.5) + dualNeq[2]*(-0.5) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(-1));

    return 1;						
}

inline int funcJ3SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=-1;
	tableau[1]=-1;
	tableau[2]=2;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=0;
	tableau[7]=0;
	tableau[8]=-1;
	
	tableau[9]=1;
	tableau[10]=2;
	tableau[11]=-1;
	tableau[12]=0;
	tableau[13]=0;
	tableau[14]=0;
	tableau[15]=-1;
	tableau[16]=0;
	tableau[17]=-1;
	
	tableau[0]=-1;
	tableau[18]=-0.5;
	tableau[19]=-0.5;
	tableau[20]=0;
	tableau[21]=0;
	tableau[22]=0;
	tableau[23]=0;
	tableau[24]=-1;
	tableau[25]=-2;
	return 1;
}


inline int funcJ3LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){  //(l2) - (Q  A^T)(y) = (c)
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
      
      matrixA[9]=0;
      matrixA[10]=0;
      matrixA[11]=0;
      
      matrixA[12]=0;
      matrixA[13]=0;
      matrixA[14]=0;
      
      matrixA[15]=-1;
      matrixA[16]=0;
      matrixA[17]=0;
      
      matrixA[18]=0;
      matrixA[19]=-1;
      matrixA[20]=0;
      
      matrixA[21]=0;
      matrixA[22]=0;
      matrixA[23]=-1;
      
      
      matrixCB[0]=1;
      matrixCB[1]=1;
      matrixCB[2]=2;
      matrixCB[3]=1;
      matrixCB[4]=-2*x[0]+1;
      matrixCB[5]=-2*x[1]+1;
      matrixCB[6]=0;
      matrixCB[7]=0;
      matrixCB[8]=0;
      matrixCB[9]=0;
      matrixCB[10]=0;

      
      return 1;
}

//const double boundJ3[10]={0,1e2,0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundJ3[10]={0,0.75,0,1.5,0,1e3,0,1.5,0,2};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J4 - solucao (10,10)*/

inline double funcJ4UP(double x[], double y[]){  //F(x,y) MAX
  return - (2*x[0]-x[1]-0.5*y[0]);
}

inline double funcJ4LW(double x[], double y[]){  //f(x,y) MAX
  return -(-x[0]-x[1]+4*y[0]-y[1]);
}

inline int funcJ4CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ4CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcJ4CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ4CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+y[1]-2.5);
  constraintValuesListReturn[1]=-(-x[0]+3*x[1]-y[1]+2);
  constraintValuesListReturn[2]=-(-x[0]-x[1]+2);
  constraintValuesListReturn[3]=-x[0];
  constraintValuesListReturn[4]=-x[1];
  constraintValuesListReturn[5]=-y[0];
  constraintValuesListReturn[6]=-y[1];
  
  return 1;
}


inline int funcJ4CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) MAX
    constraintValuesListReturn[0]=(4) - (dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1) + dualNeq[6]*(0));
    constraintValuesListReturn[1]=(-1) - (dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1));

    return 1;						
}

inline int funcJ4SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=1;
	tableau[1]=0;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=0;
	tableau[7]=4;
	
	tableau[8]=-1;
	tableau[9]=1;
	tableau[10]=0;
	tableau[11]=0;
	tableau[12]=0;
	tableau[13]=0;
	tableau[14]=-1;
	tableau[15]=-1;
	
	return 1;
}

//const double boundJ4[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundJ4[8]={0,2,0,2,0,5.5,0,8};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J5 - solucao (10,10)*/

inline double funcJ5UP(double x[], double y[]){  //F(x,y) 
  return 2*x[0]+2*x[1]-3*y[0]-3*y[1]-60;
}

inline double funcJ5LW(double x[], double y[]){  //f(x,y) 
  return (y[0]-x[0]+20)*(y[0]-x[0]+20) + (y[1]-x[1]+20)*(y[1]-x[1]+20);
}

inline int funcJ5CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ5CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+x[1]+y[0]-2*y[1]-40;

  return 1;
}

inline int funcJ5CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ5CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
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


inline int funcJ5CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) 
    constraintValuesListReturn[0]=(2*(y[0]-x[0]+20)) + (dualNeq[0]*(2) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(1) + dualNeq[8]*(0) + dualNeq[9]*(0));
    constraintValuesListReturn[1]=(2*(y[1]-x[1]+20)) + (dualNeq[0]*(0) + dualNeq[1]*(2) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(0) + dualNeq[7]*(0) + dualNeq[8]*(-1) + dualNeq[9]*(1));

    return 1;						
}

inline int funcJ5SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
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
inline int funcJ5LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
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


const double boundJ5[8]={0,50,0,50,-10,20,-10,20};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J6 - solucao (10,10)*/

inline double funcJ6UP(double x[], double y[]){  //F(x,y) 
  return (x[0]-5)*(x[0]-5)+(2*y[0]+1)*(2*y[0]+1);
}

inline double funcJ6LW(double x[], double y[]){  //f(x,y) 
  return (y[0]-1)*(y[0]-1) - 1.5*x[0]*y[0];
}

inline int funcJ6CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ6CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  return 1;
}

inline int funcJ6CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}


inline int funcJ6CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(3*x[0]-y[0]-3);
  constraintValuesListReturn[1]=-(-x[0]+0.5*y[0]+4);
  constraintValuesListReturn[2]=-(-x[0]-y[0]+7);
  constraintValuesListReturn[3]=-x[0];
  constraintValuesListReturn[4]=-y[0];
  
  return 1;
}


inline int funcJ6CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) 
    constraintValuesListReturn[0]=(2*(y[0]-1) -1.5*x[0]) + (dualNeq[0]*(1) + dualNeq[1]*(-0.5) + dualNeq[2]*(1) + dualNeq[3]*(0) + dualNeq[4]*(-1));

    return 1;						
}

inline int funcJ6SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=1;
	tableau[1]=-0.5;
	tableau[2]=1;
	tableau[3]=0;
	tableau[4]=-1;
	tableau[5]=-(2*(y[0]-1) -1.5*x[0]);
	
	return 1;
}


inline int funcJ6LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)
      
      matrixA[0]=1;
      
      matrixA[1]=-0.5;
      
      matrixA[2]=1;
      
      matrixA[3]=0;
      
      matrixA[4]=-1;
      
      
      matrixCB[0]=-1.5*x[0]-2;  
      matrixCB[1]=3*x[0]-3;
      matrixCB[2]=-x[0]+4;
      matrixCB[3]=-x[0]+7;
      matrixCB[4]=0;
      matrixCB[5]=0;
      
      return 1;
}



//const double boundJ6[8]={0,1e2,0,1e2};  //Bounds x, y
const double boundJ6[8]={0,5,0,4.5};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J7 - solucao (10,10)*/

inline double funcJ7UP(double x[], double y[]){  //F(x,y) 
  return -x[0]*x[0]-3*x[1]-4*y[0]+y[1]*y[1];
}

inline double funcJ7LW(double x[], double y[]){  //f(x,y) 
  return 2*x[0]*x[0]+y[0]*y[0]-5*y[1];
}

inline int funcJ7CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ7CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]*x[0]+2*x[1]-4;
  
  return 1;
}

inline int funcJ7CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ7CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(x[0]*x[0]-2*x[0]+x[1]*x[1]-2*y[0]+y[1]+3);
  constraintValuesListReturn[1]=-(x[1]+3*y[0]-4*y[1]-4);
  constraintValuesListReturn[2]=-x[0];
  constraintValuesListReturn[3]=-x[1];
  constraintValuesListReturn[4]=-y[0];
  constraintValuesListReturn[5]=-y[1];
  
  return 1;
}


inline int funcJ7CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y) 
    constraintValuesListReturn[0]=(2*(y[0])) + (dualNeq[0]*(2) + dualNeq[1]*(-3) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=(-5) + (dualNeq[0]*(-1) + dualNeq[1]*(4) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(0) + dualNeq[5]*(-1));

    return 1;						
}

inline int funcJ7SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=2;
	tableau[1]=-3;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=-1;
	tableau[5]=0;
	tableau[6]=-(2*(y[0]));
	
	tableau[7]=-1;
	tableau[8]=4;
	tableau[9]=0;
	tableau[10]=0;
	tableau[11]=0;
	tableau[12]=-1;
	tableau[13]=5;
	
	return 1;
}




inline int funcJ7LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=2;							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      
      matrixA[0]=2;
      matrixA[1]=-1;
      
      matrixA[2]=-3;
      matrixA[3]=4;
      
      matrixA[4]=0;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=0;

      matrixA[8]=-1;
      matrixA[9]=0;

      matrixA[10]=0;
      matrixA[11]=-1;

      
      matrixCB[0]=0;  
      matrixCB[1]=-5;
      matrixCB[2]=x[0]*x[0]-2*x[0]+x[1]*x[1]+3;
      matrixCB[3]=x[1]-4;
      matrixCB[4]=0;
      matrixCB[5]=0;
      matrixCB[6]=0;
      matrixCB[7]=0;
      
      return 1;
}


//const double boundJ7[8]={0,1e2,0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundJ7[8]={0,2,0,2,0,6,0,4};  //Bounds x, y Aproximadamente


/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J8 - solucao (10,10)*/

inline double funcJ8UP(double x[], double y[]){  //F(x,y) MAX
  return - (x[0]+3*y[0]);
}

inline double funcJ8LW(double x[], double y[]){  //f(x,y)   MAX
  return - (x[0] -3*y[0]);
}

inline int funcJ8CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ8CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ8CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ8CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]-2*y[0]+10;
  constraintValuesListReturn[1]=x[0]-2*y[0]-6;
  constraintValuesListReturn[2]=2*x[0]-y[0]-21;
  constraintValuesListReturn[3]=x[0]+2*y[0]-38;
  constraintValuesListReturn[4]=-x[0]+2*y[0]-18;
  constraintValuesListReturn[5]=-x[0];
  constraintValuesListReturn[6]=-y[0];
  
  return 1;
}


inline int funcJ8CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  MAX
    constraintValuesListReturn[0]=(-3) - (dualNeq[0]*(-2) + dualNeq[1]*(-2) + dualNeq[2]*(-1) + dualNeq[3]*(2) + dualNeq[4]*(2) + dualNeq[5]*(0) + dualNeq[6]*(-1));

    return 1;						
}

inline int funcJ8SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
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

//const double boundJ8[4]={0,1e2,0,1e2};  //Bounds x, y
const double boundJ8[4]={0,16,0,14};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J9 - solucao (10,10)*/

inline double funcJ9UP(double x[], double y[]){  //F(x,y) 
  return (x[0]-1)*(x[0]-1) + 2*y[0] - 2*x[0];
}

inline double funcJ9LW(double x[], double y[]){  //f(x,y)  
  return (2*y[0]-4)*(2*y[0]-4) + (2*y[1]-1)*(2*y[1]-1) + x[0]*y[0];
}

inline int funcJ9CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ9CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ9CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ9CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=4*x[0]+5*y[0]+4*y[1]-12;
  constraintValuesListReturn[1]=-4*x[0]-5*y[0]+4*y[1]+4;
  constraintValuesListReturn[2]=4*x[0]-4*y[0]+5*y[1]-4;
  constraintValuesListReturn[3]=-4*x[0]+4*y[0]+5*y[1]-4;
  constraintValuesListReturn[4]=-x[0];
  constraintValuesListReturn[5]=-y[0];
  constraintValuesListReturn[6]=-y[1];
  
  return 1;
}


inline int funcJ9CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*2*(2*y[0]-4) + x[0]) + (dualNeq[0]*(5) + dualNeq[1]*(-5) + dualNeq[2]*(-4) + dualNeq[3]*(4) + dualNeq[4]*(0) + dualNeq[5]*(-1) + dualNeq[6]*(0));
    constraintValuesListReturn[1]=(2*2*(2*y[1]-1)) + (dualNeq[0]*(4) + dualNeq[1]*(4) + dualNeq[2]*(5) + dualNeq[3]*(5) + dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1));

    return 1;						
}

inline int funcJ9SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=5;
	tableau[1]=-5;
	tableau[2]=-4;
	tableau[3]=4;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=0;
	tableau[7]=-(2*2*(2*y[0]-4) + x[0]);
	
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


inline int funcJ9LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=8;							//(f)  - (-A  0 )(l1)= (-b)
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

      
      matrixCB[0]=-16 + x[0];  
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

//const double boundJ9[6]={0,1e2,0,1e2,0,1e2};  //Bounds x, y
const double boundJ9[6]={0,2,0,2,0,1};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/

/* Função J10 - solucao (10,10)*/

inline double funcJ10UP(double x[], double y[]){  //F(x,y) 
  return x[0]*x[0] - 2*x[0] + x[1]*x[1]  - 2*x[1] + y[0]*y[0] + y[1]*y[1];
}

inline double funcJ10LW(double x[], double y[]){  //f(x,y)  
  return (y[0]-x[0])*(y[0]-x[0]) + (y[1]-x[1])*(y[1]-x[1]);
}

inline int funcJ10CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ10CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ10CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ10CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=-x[1];
  constraintValuesListReturn[2]=-y[0]+0.5;
  constraintValuesListReturn[3]=y[0]-1.5;
  constraintValuesListReturn[4]=-y[1]+0.5;
  constraintValuesListReturn[5]=y[1]-1.5;
  
  return 1;
}


inline int funcJ10CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(y[0]-x[0])) + (dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1) + dualNeq[3]*(1) + dualNeq[4]*(0) + dualNeq[5]*(0));
    constraintValuesListReturn[1]=(2*(y[1]-x[1])) + (dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(0) + dualNeq[3]*(0) + dualNeq[4]*(-1) + dualNeq[5]*(1));

    return 1;						
}

inline int funcJ10SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=0;
	tableau[1]=0;
	tableau[2]=-1;
	tableau[3]=1;
	tableau[4]=0;
	tableau[5]=0;
	tableau[6]=-(2*(y[0]-x[0]));
	
	tableau[7]=0;
	tableau[8]=0;
	tableau[9]=0;
	tableau[10]=0;
	tableau[11]=-1;
	tableau[12]=1;
	tableau[13]=-(2*(y[1]-x[1]));
	
	return 1;
}

const double boundJ10[8]={0,1e2,0,1e2,0.5,1.5,0.5,1.5};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J11 - solucao (10,10)*/

inline double funcJ11UP(double x[], double y[]){  //F(x,y) 
  return 16*x[0]*x[0] + 9*y[0]*y[0];
}

inline double funcJ11LW(double x[], double y[]){  //f(x,y)  
  return (x[0] + y[0] - 20)*(x[0] + y[0] - 20)*(x[0] + y[0] - 20)*(x[0] + y[0] - 20);
}

inline int funcJ11CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ11CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-4*x[0]+y[0];
  
  return 1;
}

inline int funcJ11CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ11CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=4*x[0]+y[0]-50;
  constraintValuesListReturn[1]=-x[0];
  constraintValuesListReturn[2]=-y[0];
  
  return 1;
}


inline int funcJ11CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(4*(x[0] + y[0] - 20)*(x[0] + y[0] - 20)*(x[0] + y[0] - 20)) + (dualNeq[0]*(1) + dualNeq[1]*(0) + dualNeq[2]*(-1));

    return 1;						
}

inline int funcJ11SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=1;
	tableau[1]=0;
	tableau[2]=-1;
	tableau[3]=-(4*(x[0] + y[0] - 20)*(x[0] + y[0] - 20)*(x[0] + y[0] - 20));
	
	return 1;
}

//const double boundJ11[4]={0,1e2,0,1e2};  //Bounds x, y
const double boundJ11[4]={0,12.5,0,25};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J12 - solucao (10,10)*/

inline double funcJ12UP(double x[], double y[]){  //F(x,y) 
  return x[0]-4*y[0];
}

inline double funcJ12LW(double x[], double y[]){  //f(x,y)  
  return y[0];
}

inline int funcJ12CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ12CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ12CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ12CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]-y[0]+3;
  constraintValuesListReturn[1]=-2*x[0]+y[0];
  constraintValuesListReturn[2]=2*x[0]+y[0]-12;
  constraintValuesListReturn[3]=-(-3*x[0]+2*y[0]+4);
  constraintValuesListReturn[4]=-x[0];
  constraintValuesListReturn[5]=-y[0];
  return 1;
}


inline int funcJ12CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(1) + (dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(1) + dualNeq[3]*(-2) + dualNeq[4]*(0) + dualNeq[5]*(-1));

    return 1;						
}

inline int funcJ12SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=-1;
	tableau[1]=1;
	tableau[2]=1;
	tableau[3]=-2;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=-1;
	
	return 1;
}

//const double boundJ12[4]={0,1e2,0,1e2};  //Bounds x, y
const double boundJ12[4]={0,4,0,6};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J13 - solucao (10,10)*/

inline double funcJ13UP(double x[], double y[]){  //F(x,y) 
  return x[0] + y[0];
}

inline double funcJ13LW(double x[], double y[]){  //f(x,y)  
  return -5*x[0] - y[0];
}

inline int funcJ13CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ13CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ13CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ13CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0]-0.5*y[0]+2;
  constraintValuesListReturn[1]=-0.25*x[0]+y[0]-2;
  constraintValuesListReturn[2]=x[0]+0.5*y[0]-8;
  constraintValuesListReturn[3]=x[0]-2*y[0]-4;
  constraintValuesListReturn[4]=-x[0];
  constraintValuesListReturn[5]=-y[0];
  
  return 1;
}


inline int funcJ13CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(-5) + (dualNeq[0]*(-0.5) + dualNeq[1]*(1) + dualNeq[2]*(0.5) + dualNeq[3]*(-2) + dualNeq[4]*(0) + dualNeq[5]*(-1));

    return 1;						
}

inline int funcJ13SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=-0.5;
	tableau[1]=2;
	tableau[2]=0.5;
	tableau[3]=-2;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=5;
	
	return 1;
}

//const double boundJ13[4]={0,1e2,0,1e2};  //Bounds x, y
const double boundJ13[4]={0,7.2,0,3.6};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J14 - solucao (10,10)*/

inline double funcJ14UP(double x[], double y[]){  //F(x,y) 
  return (x[0]-1)*(x[0]-1) + (y[0]-1)*(y[0]-1);
}

inline double funcJ14LW(double x[], double y[]){  //f(x,y)  
  return 0.5*y[0]*y[0] + 500*y[0] - 50*x[0]*y[0];
}

inline int funcJ14CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ14CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ14CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ14CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=-y[0];
  
  return 1;
}


inline int funcJ14CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(0.5*2*y[0]+500-50*x[0]) + (dualNeq[0]*(0) + dualNeq[1]*(-1));

    return 1;						
}

inline int funcJ14SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=0;
	tableau[1]=-1;
	tableau[2]=-(0.5*2*y[0]+500-50*x[0]);
	
	return 1;
}

inline int funcJ14LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=1;							//(f)  - (-A  0 )(l1)= (-b)
      
      matrixA[0]=0;
      
      matrixA[1]=-1;
     
      matrixCB[0]=-50*x[0]+500;  
      matrixCB[1]=0;
      matrixCB[2]=0;
     
      
      return 1;
}

const double boundJ14[4]={0,1e2,0,1e2};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J15 - solucao (10,10)*/

inline double funcJ15UP(double x[], double y[]){  //F(x,y) MAX
  return -(100*x[0]+1000*y[0]);
}

inline double funcJ15LW(double x[], double y[]){  //f(x,y)  MAX
  return -(y[0]+y[1]);
}

inline int funcJ15CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ15CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ15CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ15CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]+y[0]-y[1]-1;
  constraintValuesListReturn[1]=y[0]+y[1]-1;
  constraintValuesListReturn[2]=-x[0];
  constraintValuesListReturn[3]=x[0]-1;
  constraintValuesListReturn[4]=-y[0];
  constraintValuesListReturn[5]=y[0]-1;
  constraintValuesListReturn[6]=-y[1];
  constraintValuesListReturn[7]=y[1]-1;  
  return 1;
}


inline int funcJ15CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  MAX
    constraintValuesListReturn[0]=(1) - (dualNeq[0]*(1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0) +dualNeq[4]*(-1) + dualNeq[5]*(1) + dualNeq[6]*(0) + dualNeq[7]*(0));
    constraintValuesListReturn[1]=(1) - (dualNeq[0]*(-1) + dualNeq[1]*(1) + dualNeq[2]*(0) + dualNeq[3]*(0) +dualNeq[4]*(0) + dualNeq[5]*(0) + dualNeq[6]*(-1) + dualNeq[7]*(1));
    
    return 1;						
}

inline int funcJ15SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=1;
	tableau[1]=1;
	tableau[2]=0;
	tableau[3]=0;
	tableau[4]=-1;
	tableau[5]=1;
	tableau[6]=0;
	tableau[7]=0;
	tableau[8]=1;
	
	tableau[9]=-1;
	tableau[10]=1;
	tableau[11]=0;
	tableau[12]=0;
	tableau[13]=0;
	tableau[14]=0;
	tableau[15]=-1;
	tableau[16]=1;
	tableau[17]=1;
	
	
	return 1;
}


inline int funcJ15LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							//(f)  - (-A  0 )(l1)= (-b)
      matrixQ[1]=0;
      matrixQ[2]=0;
      matrixQ[3]=0;
      
      matrixA[0]=1;
      matrixA[1]=-1;
      
      matrixA[2]=1;
      matrixA[3]=1;
      
      matrixA[4]=0;
      matrixA[5]=0;
      
      matrixA[6]=0;
      matrixA[7]=0;
      
      matrixA[8]=-1;
      matrixA[9]=0;
      
      matrixA[10]=1;
      matrixA[11]=0;
      
      matrixA[12]=0;
      matrixA[13]=-1;
      
      matrixA[14]=0;
      matrixA[15]=1;
      
     
      matrixCB[0]=-1;  
      matrixCB[1]=-1;
      matrixCB[2]=-x[0]+1;
      matrixCB[3]=1;
      matrixCB[4]=0;
      matrixCB[5]=1;
      matrixCB[6]=0;
      matrixCB[7]=1;
      matrixCB[8]=0;
      matrixCB[9]=1;
     
      
      return 1;
}


const double boundJ15[6]={0,1,0,1,0,1};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J16 - solucao (10,10)*/

inline double funcJ16UP(double x[], double y[]){  //F(x,y) 
  return (x[0]-3)*(x[0]-3) + (y[0]-2)*(y[0]-2);
}

inline double funcJ16LW(double x[], double y[]){  //f(x,y)  
  return (y[0]-5)*(y[0]-5);
}

inline int funcJ16CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ16CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ16CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ16CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+1);
  constraintValuesListReturn[1]=-(-x[0]+2*y[0]-2);
  constraintValuesListReturn[2]=-(-x[0]-2*y[0]+14);
  constraintValuesListReturn[3]=-x[0];
  constraintValuesListReturn[4]=x[0]-8;
  constraintValuesListReturn[5]=-y[0];
  
  return 1;
}


inline int funcJ16CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(y[0]-5)) + (dualNeq[0]*(1) + dualNeq[1]*(-2) + dualNeq[2]*(2) + dualNeq[3]*(0) +dualNeq[4]*(0) + dualNeq[5]*(-1));
    
    return 1;						
}

inline int funcJ16SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=1;
	tableau[1]=-2;
	tableau[2]=2;
	tableau[3]=0;
	tableau[4]=0;
	tableau[5]=-1;
	tableau[6]=-(2*(y[0]-5));	
	
	return 1;
}

//const double boundJ16[4]={0,8,0,1e2};  //Bounds x, y
const double boundJ16[4]={0,8,0,5.8};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/



/* Função J17 - solucao (10,10)*/

inline double funcJ17UP(double x[], double y[]){  //F(x,y) 
  return (x[0]-3)*(x[0]-3) + (y[0]-2)*(y[0]-2);
}

inline double funcJ17LW(double x[], double y[]){  //f(x,y)  
  return (y[0]-5)*(y[0]-5);
}

inline int funcJ17CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ17CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  constraintValuesListReturn[0]=-(2*x[0]-y[0]+1);
  constraintValuesListReturn[1]=-(-x[0]+2*y[0]-2);
  constraintValuesListReturn[2]=-(-x[0]-2*y[0]+14);
  return 1;
}

inline int funcJ17CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ17CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=-x[0];
  constraintValuesListReturn[1]=x[0]-8;
  constraintValuesListReturn[2]=-y[0];
  
  return 1;
}


inline int funcJ17CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  
    constraintValuesListReturn[0]=(2*(y[0]-5)) + (dualNeq[0]*(0) + dualNeq[1]*(0) + dualNeq[2]*(-1));
    
    return 1;						
}

inline int funcJ17SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = - grad(f(x,y))  
	tableau[0]=0;
	tableau[1]=0;
	tableau[2]=-1;
	tableau[3]=-(2*(y[0]-5));	
	
	return 1;
}

//const double boundJ17[4]={0,8,0,1e2};  //Bounds x, y
const double boundJ17[4]={0,8,0,5.8};  //Bounds x, y


/* --------------------------------------------------------------------------------------------------------------------------------------*/


/* Função J18 - solucao (10,10)*/

inline double funcJ18UP(double x[], double y[]){  //F(x,y) MAX
  return -(-2*x[0]+11*y[0]);
}

inline double funcJ18LW(double x[], double y[]){  //f(x,y)  MAX
  return -(-x[0]-3*y[0]);
}

inline int funcJ18CTREQUP(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
    return 1;
}

inline int funcJ18CTRNEQUP(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0  
  return 1;
}

inline int funcJ18CTREQLW(double x[], double y[], double constraintValuesListReturn[]){ //g(x,y)=0
  return 1;
}


inline int funcJ18CTRNEQLW(double x[], double y[], double constraintValuesListReturn[]){  //g(x,y)<=0
  constraintValuesListReturn[0]=x[0]-2*y[0]-4;
  constraintValuesListReturn[1]=2*x[0]-y[0]-24;
  constraintValuesListReturn[2]=3*x[0]+4*y[0]-96;
  constraintValuesListReturn[3]=x[0]+7*y[0]-126;//constraintValuesListReturn[3]=x[0]+4*y[0]-126;
  constraintValuesListReturn[4]=-4*x[0]+5*y[0]-65;
  constraintValuesListReturn[5]=-(x[0]+4*y[0]-8);
  constraintValuesListReturn[6]=-x[0];
  constraintValuesListReturn[7]=-y[0];
 
  return 1;
}


inline int funcJ18CTKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){ //grad Lagrangeano(x,y)  MAX
    constraintValuesListReturn[0]=(-3) - (dualNeq[0]*(-2) + dualNeq[1]*(-1) + dualNeq[2]*(4) + dualNeq[3]*(4) + dualNeq[4]*(5) + dualNeq[5]*(-4) + dualNeq[6]*(0) + dualNeq[7]*(-1));
    
    return 1;						
}

inline int funcJ18SimplexTableauKKT(double x[], double y[], double tableau[]){  //grad(h(x,y)) \lambda = grad(f(x,y))  MAX
	tableau[0]=-2;
	tableau[1]=-1;
	tableau[2]=4;
	tableau[3]=4;
	tableau[4]=5;
	tableau[5]=-4;
	tableau[6]=0;
	tableau[7]=-1;
	tableau[8]=-3;
	
	return 1;
}



inline int funcJ18LemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){//(l2) - (Q  A^T)(y) = (c)
      matrixQ[0]=0;							//(f)  - (-A  0 )(l1)= (-b)

      
      matrixA[0]=-2;
      
      matrixA[1]=-1;
      
      matrixA[2]=4;
      
      matrixA[3]=7;//matrixA[3]=4;
      
      matrixA[4]=5;
      
      matrixA[5]=-4;
      
      matrixA[6]=0;
      
      matrixA[7]=-1;
      
      
     
      matrixCB[0]=3;  
      matrixCB[1]=-x[0]+4;
      matrixCB[2]=-2*x[0]+24;
      matrixCB[3]=-3*x[0]+96;
      matrixCB[4]=-x[0]+126;
      matrixCB[5]=4*x[0]+65;
      matrixCB[6]=x[0]-8;
      matrixCB[7]=0;
      matrixCB[8]=0;
     
      
      return 1;
}

//const double boundJ18[4]={-1e2,0,-1e2,0};  //Bounds x, y
const double boundJ18[4]={0,18,0,18};  //Bounds x, y

/* --------------------------------------------------------------------------------------------------------------------------------------*/

const FunctionJAQ listFunctionJaq[DEFINEfunctionListSizeJaq]={{1,1,0,3,0,3,boundJ1,funcJ1UP,funcJ1LW,funcJ1CTREQUP,funcJ1CTRNEQUP, funcJ1CTREQLW,funcJ1CTRNEQLW,funcJ1CTKKT,funcJ1SimplexTableauKKT,funcJ1LemkeMatrix,"funcJ1"},
					  {2,2,0,3,0,4,boundJ2,funcJ2UP,funcJ2LW,funcJ2CTREQUP,funcJ2CTRNEQUP, funcJ2CTREQLW,funcJ2CTRNEQLW,funcJ2CTKKT,funcJ2SimplexTableauKKT,funcJ2LemkeMatrix,"funcJ2"},
					  {2,3,0,0,0,8,boundJ3,funcJ3UP,funcJ3LW,funcJ3CTREQUP,funcJ3CTRNEQUP, funcJ3CTREQLW,funcJ3CTRNEQLW,funcJ3CTKKT,funcJ3SimplexTableauKKT,funcJ3LemkeMatrix,"funcJ3"},
					  {2,2,0,0,0,7,boundJ4,funcJ4UP,funcJ4LW,funcJ4CTREQUP,funcJ4CTRNEQUP, funcJ4CTREQLW,funcJ4CTRNEQLW,funcJ4CTKKT,funcJ4SimplexTableauKKT,NULL,"funcJ4"},
					  {2,2,0,1,0,10,boundJ5,funcJ5UP,funcJ5LW,funcJ5CTREQUP,funcJ5CTRNEQUP, funcJ5CTREQLW,funcJ5CTRNEQLW,funcJ5CTKKT,funcJ5SimplexTableauKKT,funcJ5LemkeMatrix,"funcJ5"},
					  {1,1,0,0,0,5,boundJ6,funcJ6UP,funcJ6LW,funcJ6CTREQUP,funcJ6CTRNEQUP, funcJ6CTREQLW,funcJ6CTRNEQLW,funcJ6CTKKT,funcJ6SimplexTableauKKT,funcJ6LemkeMatrix,"funcJ6"},
					  {2,2,0,1,0,6,boundJ7,funcJ7UP,funcJ7LW,funcJ7CTREQUP,funcJ7CTRNEQUP, funcJ7CTREQLW,funcJ7CTRNEQLW,funcJ7CTKKT,funcJ7SimplexTableauKKT,funcJ7LemkeMatrix,"funcJ7"},
					  {1,1,0,0,0,7,boundJ8,funcJ8UP,funcJ8LW,funcJ8CTREQUP,funcJ8CTRNEQUP, funcJ8CTREQLW,funcJ8CTRNEQLW,funcJ8CTKKT,funcJ8SimplexTableauKKT,NULL,"funcJ8"},
					  {1,2,0,0,0,7,boundJ9,funcJ9UP,funcJ9LW,funcJ9CTREQUP,funcJ9CTRNEQUP, funcJ9CTREQLW,funcJ9CTRNEQLW,funcJ9CTKKT,funcJ9SimplexTableauKKT,funcJ9LemkeMatrix,"funcJ9"},
					  {2,2,0,0,0,6,boundJ10,funcJ10UP,funcJ10LW,funcJ10CTREQUP,funcJ10CTRNEQUP, funcJ10CTREQLW,funcJ10CTRNEQLW,funcJ10CTKKT,funcJ10SimplexTableauKKT,NULL,"funcJ10"},
					  {1,1,0,1,0,3,boundJ11,funcJ11UP,funcJ11LW,funcJ11CTREQUP,funcJ11CTRNEQUP, funcJ11CTREQLW,funcJ11CTRNEQLW,funcJ11CTKKT,funcJ11SimplexTableauKKT,NULL,"funcJ11"},
					  {1,1,0,0,0,6,boundJ12,funcJ12UP,funcJ12LW,funcJ12CTREQUP,funcJ12CTRNEQUP, funcJ12CTREQLW,funcJ12CTRNEQLW,funcJ12CTKKT,funcJ12SimplexTableauKKT,NULL,"funcJ12"},
					  {1,1,0,0,0,6,boundJ13,funcJ13UP,funcJ13LW,funcJ13CTREQUP,funcJ13CTRNEQUP, funcJ13CTREQLW,funcJ13CTRNEQLW,funcJ13CTKKT,funcJ13SimplexTableauKKT,NULL,"funcJ13"},
					  {1,1,0,0,0,2,boundJ14,funcJ14UP,funcJ14LW,funcJ14CTREQUP,funcJ14CTRNEQUP, funcJ14CTREQLW,funcJ14CTRNEQLW,funcJ14CTKKT,funcJ14SimplexTableauKKT,funcJ14LemkeMatrix,"funcJ14"},
					  {1,2,0,0,0,8,boundJ15,funcJ15UP,funcJ15LW,funcJ15CTREQUP,funcJ15CTRNEQUP, funcJ15CTREQLW,funcJ15CTRNEQLW,funcJ15CTKKT,funcJ15SimplexTableauKKT,funcJ15LemkeMatrix,"funcJ15"},
					  {1,1,0,0,0,6,boundJ16,funcJ16UP,funcJ16LW,funcJ16CTREQUP,funcJ16CTRNEQUP, funcJ16CTREQLW,funcJ16CTRNEQLW,funcJ16CTKKT,funcJ16SimplexTableauKKT,NULL,"funcJ16"},
					  {1,1,0,3,0,3,boundJ17,funcJ17UP,funcJ17LW,funcJ17CTREQUP,funcJ17CTRNEQUP, funcJ17CTREQLW,funcJ17CTRNEQLW,funcJ17CTKKT,funcJ17SimplexTableauKKT,NULL,"funcJ17"},
					  {1,1,0,0,0,8,boundJ18,funcJ18UP,funcJ18LW,funcJ18CTREQUP,funcJ18CTRNEQUP, funcJ18CTREQLW,funcJ18CTRNEQLW,funcJ18CTKKT,funcJ18SimplexTableauKKT,funcJ18LemkeMatrix,"funcJ18"}

};



#endif
