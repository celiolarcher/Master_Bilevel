#ifndef SOLUTION_INCLUDED
#define SOLUTION_INCLUDED

#include "InputFunction.h"

class Solution{
    public: double *vet;
    public: int sizeVet;
  
    public: double calcScore(InputFunction *function);
    public: int initValue(double initVet[], InputFunction *function);
    public: int initRandom(InputFunction *function);
    
    Solution(InputFunction *function){
        sizeVet=function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();
        vet=new double[sizeVet];
    }
    
    ~Solution(){
        delete vet;
    }
};


#endif
