#ifndef INPUTFUNCTION_INCLUDED
#define INPUTFUNCTION_INCLUDED
#include "ListFunction.h"

class InputFunction{
    private: double (*objFuncUPLevel)(double x[], double y[]);
    private: double (*objFuncLWLevel)(double x[], double y[]);
    private: int (*constrEqUP)(double x[], double y[], double constraintValuesListReturn[]);
    private: int (*constrNeqUP)(double x[], double y[], double constraintValuesListReturnt[]);
    private: int (*constrEqLW)(double x[], double y[], double constraintValuesListReturn[]);
    private: int (*constrNeqLW)(double x[], double y[], double constraintValuesListReturnt[]);
    private: int (*constrKKT)(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]);
    private: int dimensionUP;
	 int dimensionLW;
	 int countEqConstraintUP;
	 int countNeqConstraintUP;
	 int countEqConstraintLW;
	 int countNeqConstraintLW;	
    public:  const double *bounds;
    
    public: InputFunction(char *funcName);
    public: double getUPLevelFunction(double x[], double y[]);
    public: double getLWLevelFunction(double x[], double y[]);
    public: int constraintsValueEQUP(double x[], double y[], double constraintValuesListReturn[]);  //Válida se =0
    public: int constraintsValueNEQUP(double x[], double y[], double constraintValuesListReturn[]); //Válida se <=0
    public: int constraintsValueEQLW(double x[], double y[], double constraintValuesListReturn[]);  //Válida se =0
    public: int constraintsValueNEQLW(double x[], double y[], double constraintValuesListReturn[]); //Válida se <=0
    public: int constraintsValueKKT(double x[], double y[], double dualEq[],  double  dualNeq[], double constraintValuesListReturn[]); //Válida se =0
    public: int constraintsSlackness(double  dualNeq[], double constraintNeqValueList[], double constraintValuesListReturn[]); //Válida se =0
    public: int constraintsSlackness(double x[], double y[], double  dualNeq[], double constraintValuesListReturn[]); //Válida se =0
    public: int setFunction(char *nameFunc);    
  
    
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
};


#endif