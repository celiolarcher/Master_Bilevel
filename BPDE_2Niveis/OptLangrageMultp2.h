#ifndef PENALTYFUNCTION_INCLUDED
#define PENALTYFUNCTION_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LagrangeMultpAPM: public SolutionDecoder{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
	int initInstance(InputFunction *function);
        int decodifySolution(Solution *sol);
        int updatePenalty(Solution **population, Solution **newPop, int sizePop, int sizeNewPop);
        int compareSolutions(Solution *sol1, Solution *sol2);
};

double validateSolution(Solution *sol,InputFunction *function);
int compareSolutions(Solution *sol1, Solution *sol2,InputFunction *function);
int updateAPM(Solution **population, int sizePop, InputFunction *function);
int compareAPM(Solution *sol1, Solution *sol2,InputFunction *function);

int updateAPMSmooth(Solution **population, int sizePop, InputFunction *function);
int compareAPMSmooth(Solution *sol1, Solution *sol2,InputFunction *function);
int updateAPMSmoothCHKS(Solution **population, int sizePop, InputFunction *function);

double getLambda(Solution *sol, int pos, InputFunction *function);
#endif
 
