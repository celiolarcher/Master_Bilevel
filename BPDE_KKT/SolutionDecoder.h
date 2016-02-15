#ifndef SOLUTIONDECODER_INCLUDED
#define SOLUTIONDECODER_INCLUDED

#include "Solution.h"
#include "InputFunction.h"

class SolutionDecoder{
  public:
      int solutionSize;
      int constraintsNumber;
      double *boundAttributes;
      InputFunction *function;
  public:
      virtual int initInstance(InputFunction *function)=0;
      virtual int decodifySolution(Solution *sol)=0;
      virtual int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop)=0;
      virtual int compareSolutions(Solution *sol1, Solution *sol2)=0;
}; 

#endif