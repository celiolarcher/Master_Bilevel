#ifndef SOLUTION_INCLUDED
#define SOLUTION_INCLUDED

//#define TOL_EQ_CONST 10e-3
//#define TOL_NEQ_CONST 10e-5

#include "InputFunction.h"
#include <iostream>


extern double TOL_EQ_CONST;
extern double TOL_NEQ_CONST;

class Solution{
    public: double *vectorCharacters;
    public: double *constraintValues;
    public: int sizeVec;
    public: int countConstraint;
    public: bool feasible;
    public: double upLevelFunction;
    public: double score;
    
    public: double calcScore(InputFunction *function);
    public: int initValue(double initVec[]);
    public: int initRandom(double bounds[]);
    public:	Solution *clone();
    public: double diffZeroSolution(Solution *sol);
    public: double diffSquareSolution(Solution *sol);
    public: double diffMaxSolution(Solution *sol);
    friend  std::ostream& operator<<(std::ostream &out, Solution &sol);  //Imprime a configuração da solução.
    
    Solution(int size, int constraints){
        sizeVec=size;
        countConstraint=constraints;
        vectorCharacters=new double[sizeVec];
        if(countConstraint>0)
          constraintValues=new double[countConstraint];
        else constraintValues=NULL;
    }
    
    ~Solution(){
        delete vectorCharacters;
        if(constraintValues) delete constraintValues;
    }
    
};


#endif
  
