#ifndef LAGRANGEMULTPSMOOTH_INCLUDED
#define LAGRANGEMULTPSMOOTH_INCLUDED  

#include "SolutionDecoder.h"

class LagrangeMultpSmooth: public SolutionDecoder{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
	int initInstance(InputFunction *function);
        int decodifySolution(Solution *sol);
};

#endif
 
