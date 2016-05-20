#ifndef APMPENALTY_INCLUDED
#define APMPENALTY_INCLUDED  

#include "Penalty.h"

class APMPenalty: public PenaltySolution{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
        int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop, SolutionDecoder *decoder);
        int compareSolutions(Solution *sol1, Solution *sol2, SolutionDecoder *decoder);
};

#endif
 
