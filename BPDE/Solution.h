#ifndef SOLUTION_INCLUDED
#define SOLUTION_INCLUDED

#include "InputFunction.h"
#include <iostream>

class Solution{
    public: double *vectorCharacters;
    public: int sizeVec;
  
    public: double calcScore(InputFunction *function);
    public: int initValue(double initVec[]);
    public: int initRandom(double bounds[]);
    friend  std::ostream& operator<<(std::ostream &out, Solution &sol);  //Imprime a configuração da solução.
    
    Solution(int size){
        sizeVec=size;
        vectorCharacters=new double[sizeVec];
    }
    
    ~Solution(){
        delete vectorCharacters;
    }
};


#endif
