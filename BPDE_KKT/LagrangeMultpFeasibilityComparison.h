#ifndef LAGRANGEMULTPAPM_INCLUDED
#define LAGRANGEMULTPAPM_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LagrangeMultpFeasibleLevel: public SolutionDecoder{
    public:
	int initInstance(InputFunction *function);
	int decodifySolution(Solution *sol);
	int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop);
	int compareSolutions(Solution *sol1, Solution *sol2);
};

#endif
 
