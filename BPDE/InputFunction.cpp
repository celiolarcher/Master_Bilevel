#include "InputFunction.h"
#include <string.h>
    InputFunction::InputFunction(char *funcName){
         setFunction(funcName);
    }
    
    
    int InputFunction::setFunction(char* nameFunc){
      for(int i=0;i<DEFINEfunctionListSize;i++){
	if(!strcmp(nameFunc,listFunction[i].name)){
	    objFuncUPLevel=listFunction[i].funcUP;
	    objFuncLWLevel=listFunction[i].funcLW;
	    constrEqUP=listFunction[i].funcCTREQUP;
	    constrNeqUP=listFunction[i].funcCTRNEQUP;
	    constrEqLW=listFunction[i].funcCTREQLW;
	    constrNeqLW=listFunction[i].funcCTRNEQLW;
	    constrKKT=listFunction[i].funcCTRKKT;
	    dimensionUP=listFunction[i].dimensionUP;
	    dimensionLW=listFunction[i].dimensionLW;
	    countEqConstraintUP=listFunction[i].numEqConstrUP;
	    countNeqConstraintUP=listFunction[i].numNeqConstrUP;
	    countEqConstraintLW=listFunction[i].numEqConstrLW;
	    countNeqConstraintLW=listFunction[i].numNeqConstrLW;
	    bounds=listFunction[i].boundsVar;
	    return 1;
	}
      }
      return 0;
    }
    
    double InputFunction::getUPLevelFunction(double x[], double y[]){
        return (*objFuncUPLevel)(x,y);
    }
    
    double InputFunction::getLWLevelFunction(double x[], double y[]){
        return (*objFuncLWLevel)(x,y);
    }

    int InputFunction::constraintsValueEQUP(double x[], double y[], double constraintValuesList[]){
        return (*constrEqUP)(x,y,constraintValuesList);
    }


    int InputFunction::constraintsValueNEQUP(double x[], double y[], double constraintValuesList[]){
        return (*constrNeqUP)(x,y,constraintValuesList);
    }
    
    int InputFunction::constraintsValueEQLW(double x[], double y[], double constraintValuesList[]){
        return (*constrEqLW)(x,y,constraintValuesList);
    }


    int InputFunction::constraintsValueNEQLW(double x[], double y[], double constraintValuesList[]){
        return (*constrNeqLW)(x,y,constraintValuesList);
    }
    
    int InputFunction::constraintsValueKKT(double x[], double y[], double dualEq[], double  dualNeq[], double constraintValuesListReturn[]){
        return (*constrKKT)(x,y,dualEq, dualNeq, constraintValuesListReturn);
    }
    
     int InputFunction::constraintsSlackness(double  dualNeq[], double constraintNeqValueList[], double constraintValuesListReturn[]){
       for(int i=0;i<countNeqConstraintLW;i++){
         constraintValuesListReturn[i]=dualNeq[i]*constraintNeqValueList[i];
       }
       
       return 1;
    }
    
#include <iostream>
using namespace std;
    int InputFunction::constraintsSlackness(double x[], double y[], double  dualNeq[], double constraintValuesListReturn[]){
       if(!constraintsValueNEQLW(x,y,constraintValuesListReturn)) return 0;


       for(int i=0;i<countNeqConstraintLW;i++){
          constraintValuesListReturn[i]*=dualNeq[i];
                    //cout<<constraintValuesListReturn[i]<<"\n";

       }
       
       return 1;
    }