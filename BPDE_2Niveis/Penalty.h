#ifndef PENALTY_INCLUDED
#define PENALTY_INCLUDED

#include "Solution.h"
#include "SolutionDecoder.h"

class PenaltySolution{
  public:
      virtual int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop, SolutionDecoder *decoder)=0;
      virtual int compareSolutions(Solution *sol1, Solution *sol2, SolutionDecoder *decoder)=0;
}; 

#endif
