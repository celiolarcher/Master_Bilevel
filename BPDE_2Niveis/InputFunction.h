#ifndef INPUTFUNCTION_INCLUDED
#define INPUTFUNCTION_INCLUDED
#include "ListFunction.h"
#include "ListFunctionSMD.h"
#include "ListFunctionJaq.h"
#include "ListFunctionOthers.h"

class InputFunction{
    private: double (*objFuncUPLevel)(double x[], double y[]);
    private: double (*objFuncLWLevel)(double x[], double y[]);
    private: int (*constrEqUP)(double x[], double y[], double constraintValuesListReturn[]);
    private: int (*constrNeqUP)(double x[], double y[], double constraintValuesListReturnt[]);
    private: int (*constrEqLW)(double x[], double y[], double constraintValuesListReturn[]);
    private: int (*constrNeqLW)(double x[], double y[], double constraintValuesListReturnt[]);
    private: int (*constrKKT)(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]);
    private: int (*simplexTableauKKT)(double x[], double y[], double tableau[]);
    private: int (*lemkeMatrix)(double x[], double matrixQ[], double matrixA[], double matrixCB[]);
    private: int dimensionUP;
	 int dimensionLW;
	 int countEqConstraintUP;
	 int countNeqConstraintUP;
	 int countEqConstraintLW;
	 int countNeqConstraintLW;
	 int upLevelCalls;
	 int lwSimplexCalls;
    public: const double *bounds;
    
    public: InputFunction(char *funcName);
    public: InputFunction(char *funcName,int P, int Q, int R, int S);
    public: double getUPLevelFunction(double x[], double y[]);
    public: double getLWLevelFunction(double x[], double y[]);
    public: int constraintsValueEQUP(double x[], double y[], double constraintValuesListReturn[]);  //Válida se =0
    public: int constraintsValueNEQUP(double x[], double y[], double constraintValuesListReturn[]); //Válida se <=0
    public: int constraintsValueEQLW(double x[], double y[], double constraintValuesListReturn[]);  //Válida se =0
    public: int constraintsValueNEQLW(double x[], double y[], double constraintValuesListReturn[]); //Válida se <=0
    public: int constraintsValueKKT(double x[], double y[], double dualEq[],  double  dualNeq[], double constraintValuesListReturn[]); //Válida se =0
    public: int constraintsSlackness(double  dualNeq[], double constraintNeqValueList[], double constraintValuesListReturn[]); //Válida se =0
    public: int constraintsSlackness(double x[], double y[], double  dualNeq[], double constraintValuesListReturn[]); //Válida se =0
    public: int getSimplexTableauKKT(double x[], double y[], double tableau[]);
    public: int getLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]);
    public: int setFunction(char *nameFunc);    
    public: int setFunctionSMD(char *nameFunc,int P, int Q, int R, int S);
    
    public: inline int getDimensionUP(){
        return dimensionUP;
    }
    
    public: inline int getDimensionLW(){
        return dimensionLW;
    }
    
    public: inline int getEQConstraintNumberUP(){
        return countEqConstraintUP;
    }
    
    
    public: inline int getNEQConstraintNumberUP(){
        return countNeqConstraintUP;
    }
    
    public: inline int getEQConstraintNumberLW(){
        return countEqConstraintLW;
    }
    
    
    public: inline int getNEQConstraintNumberLW(){
        return countNeqConstraintLW;
    }
    
    
    public: inline int getKKTConstraintNumber(){  //Dimensão do KKT equivale ao número de variáveis na função em baixo nivel
        return dimensionLW;
    }
    
    public: inline int getUPLevelCalls(){  //Dimensão do KKT equivale ao número de variáveis na função em baixo nivel
        return upLevelCalls;
    }
    
    public: inline int getLWLevelSimplexCalls(){  //Dimensão do KKT equivale ao número de variáveis na função em baixo nivel
        return lwSimplexCalls;
    }
};


#endif
