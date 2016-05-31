#include "InputFunction.h"
#include <string.h>
    InputFunction::InputFunction(char *funcName){
         setFunction(funcName);
         upLevelCalls=0;
         lwSimplexCalls=0;
    }
    
    InputFunction::InputFunction(char *funcName, int P, int Q, int R, int S){
         setFunctionSMD(funcName,P,Q,R,S);
         upLevelCalls=0;
         lwSimplexCalls=0;
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
	    simplexTableauKKT=listFunction[i].funcSimplexTableauKKT;
	    lemkeMatrix=listFunction[i].funcLemkeMatrix;
	    return 1;
	}
      }
      
      for(int i=0;i<DEFINEfunctionListSizeJaq;i++){
	if(!strcmp(nameFunc,listFunctionJaq[i].name)){
	    objFuncUPLevel=listFunctionJaq[i].funcUP;
	    objFuncLWLevel=listFunctionJaq[i].funcLW;
	    constrEqUP=listFunctionJaq[i].funcCTREQUP;
	    constrNeqUP=listFunctionJaq[i].funcCTRNEQUP;
	    constrEqLW=listFunctionJaq[i].funcCTREQLW;
	    constrNeqLW=listFunctionJaq[i].funcCTRNEQLW;
	    constrKKT=listFunctionJaq[i].funcCTRKKT;
	    dimensionUP=listFunctionJaq[i].dimensionUP;
	    dimensionLW=listFunctionJaq[i].dimensionLW;
	    countEqConstraintUP=listFunctionJaq[i].numEqConstrUP;
	    countNeqConstraintUP=listFunctionJaq[i].numNeqConstrUP;
	    countEqConstraintLW=listFunctionJaq[i].numEqConstrLW;
	    countNeqConstraintLW=listFunctionJaq[i].numNeqConstrLW;
	    bounds=listFunctionJaq[i].boundsVar;
	    simplexTableauKKT=listFunctionJaq[i].funcSimplexTableauKKT;
	    lemkeMatrix=listFunctionJaq[i].funcLemkeMatrix;
	    return 1;
	}
      }
      
      for(int i=0;i<DEFINEfunctionListSizeOthers;i++){
	if(!strcmp(nameFunc,listFunctionOthers[i].name)){
	    objFuncUPLevel=listFunctionOthers[i].funcUP;
	    objFuncLWLevel=listFunctionOthers[i].funcLW;
	    constrEqUP=listFunctionOthers[i].funcCTREQUP;
	    constrNeqUP=listFunctionOthers[i].funcCTRNEQUP;
	    constrEqLW=listFunctionOthers[i].funcCTREQLW;
	    constrNeqLW=listFunctionOthers[i].funcCTRNEQLW;
	    constrKKT=listFunctionOthers[i].funcCTRKKT;
	    dimensionUP=listFunctionOthers[i].dimensionUP;
	    dimensionLW=listFunctionOthers[i].dimensionLW;
	    countEqConstraintUP=listFunctionOthers[i].numEqConstrUP;
	    countNeqConstraintUP=listFunctionOthers[i].numNeqConstrUP;
	    countEqConstraintLW=listFunctionOthers[i].numEqConstrLW;
	    countNeqConstraintLW=listFunctionOthers[i].numNeqConstrLW;
	    bounds=listFunctionOthers[i].boundsVar;
	    simplexTableauKKT=listFunctionOthers[i].funcSimplexTableauKKT;
	    lemkeMatrix=listFunctionOthers[i].funcLemkeMatrix;
	    return 1;
	}
      }
      return 0;
    }
    
    
    int inputP, inputQ,inputR,inputS;

    int InputFunction::setFunctionSMD(char* nameFunc, int P, int Q, int R, int S){
      for(int i=0;i<DEFINEfunctionListSizeSMD;i++){
	if(!strcmp(nameFunc,listFunctionSMD[i].name)){
	    FunctionSMD *functionSelected;
	    
	    inputP=P,inputQ=Q,inputR=R,inputS=S;	    
	    
	    listFunctionSMD[i].setFuncSMD(&functionSelected);
	    
	    objFuncUPLevel=functionSelected->funcUP;
	    objFuncLWLevel=functionSelected->funcLW;
	    constrEqUP=functionSelected->funcCTREQUP;
	    constrNeqUP=functionSelected->funcCTRNEQUP;
	    constrEqLW=functionSelected->funcCTREQLW;
	    constrNeqLW=functionSelected->funcCTRNEQLW;
	    constrKKT=functionSelected->funcCTRKKT;
	    dimensionUP=functionSelected->dimensionUP;
	    dimensionLW=functionSelected->dimensionLW;
	    countEqConstraintUP=functionSelected->numEqConstrUP;
	    countNeqConstraintUP=functionSelected->numNeqConstrUP;
	    countEqConstraintLW=functionSelected->numEqConstrLW;
	    countNeqConstraintLW=functionSelected->numNeqConstrLW;
	    bounds=functionSelected->boundsVar;
	    simplexTableauKKT=functionSelected->funcSimplexTableauKKT;
	    lemkeMatrix=functionSelected->funcLemkeMatrix;

	    delete functionSelected;
	    
	    return 1;
	}
      }

      for(int i=0;i<DEFINEfunctionListSizeNewSMD;i++){
	if(!strcmp(nameFunc,listFunctionNewSMD[i].name)){
	    FunctionNewSMD *functionSelected;
	    
	    inputP=P,inputQ=Q,inputR=R,inputS=S;	    
	    
	    listFunctionNewSMD[i].setFuncSMD(&functionSelected);
	    
	    objFuncUPLevel=functionSelected->funcUP;
	    objFuncLWLevel=functionSelected->funcLW;
	    constrEqUP=functionSelected->funcCTREQUP;
	    constrNeqUP=functionSelected->funcCTRNEQUP;
	    constrEqLW=functionSelected->funcCTREQLW;
	    constrNeqLW=functionSelected->funcCTRNEQLW;
	    constrKKT=functionSelected->funcCTRKKT;
	    dimensionUP=functionSelected->dimensionUP;
	    dimensionLW=functionSelected->dimensionLW;
	    countEqConstraintUP=functionSelected->numEqConstrUP;
	    countNeqConstraintUP=functionSelected->numNeqConstrUP;
	    countEqConstraintLW=functionSelected->numEqConstrLW;
	    countNeqConstraintLW=functionSelected->numNeqConstrLW;
	    bounds=functionSelected->boundsVar;
	    simplexTableauKKT=functionSelected->funcSimplexTableauKKT;
	    lemkeMatrix=functionSelected->funcLemkeMatrix;

	    delete functionSelected;
	    
	    return 1;
	}
      }

      return 0;
    }
    
    double InputFunction::getUPLevelFunction(double x[], double y[]){
        upLevelCalls++;
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

    int InputFunction::getSimplexTableauKKT(double x[], double y[], double tableau[]){
        lwSimplexCalls++;
        return (*simplexTableauKKT)(x,y,tableau);
    }
    
    
    int InputFunction::getLemkeMatrix(double x[], double matrixQ[], double matrixA[], double matrixCB[]){
        lwSimplexCalls++;
        return (*lemkeMatrix)(x,matrixQ,matrixA,matrixCB);
    }
    
     int InputFunction::constraintsSlackness(double  dualNeq[], double constraintNeqValueList[], double constraintValuesListReturn[]){
       for(int i=0;i<countNeqConstraintLW;i++){
         constraintValuesListReturn[i]=dualNeq[i]*constraintNeqValueList[i];
       }
       
       return 1;
    }
    

    int InputFunction::constraintsSlackness(double x[], double y[], double  dualNeq[], double constraintValuesListReturn[]){
       if(!constraintsValueNEQLW(x,y,constraintValuesListReturn)) return 0;


       for(int i=0;i<countNeqConstraintLW;i++){
          constraintValuesListReturn[i]*=dualNeq[i];
       }
       
       return 1;
    }
