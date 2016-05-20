#ifndef SOLUTIONDECODER_INCLUDED
#define SOLUTIONDECODER_INCLUDED

#include "Solution.h"
#include "InputFunction.h"

class SolutionDecoder{
  public:
      int solutionSize;
      int constraintsNumber;
      int constraintsEQNumber;
      int constraintsNEQNumber;
      int editBegin;
      int editSize;
      double *boundAttributes;
      InputFunction *function;
  public:
      virtual int initInstance(InputFunction *function)=0;
      virtual int decodifySolution(Solution *sol)=0;
}; 

#endif
