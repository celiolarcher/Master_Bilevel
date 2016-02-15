#ifndef LAGRANGEMULTPAPMSMOOTH_INCLUDED
#define LAGRANGEMULTPAPMSMOOTH_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LagrangeMultpAPMSmooth: public SolutionDecoder{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
	int initInstance(InputFunction *function);
        int decodifySolution(Solution *sol);
        int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop);
        int compareSolutions(Solution *sol1, Solution *sol2);
};

#endif
 
