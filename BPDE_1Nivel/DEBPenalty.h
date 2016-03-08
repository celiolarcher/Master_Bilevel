#ifndef DEBPENALTY_INCLUDED
#define DEBPENALTY_INCLUDED  

#include "Penalty.h"

class DEBPenalty: public PenaltySolution{
    public:
        int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop, SolutionDecoder *decoder);
        int compareSolutions(Solution *sol1, Solution *sol2, SolutionDecoder *decoder);
};

#endif
 
