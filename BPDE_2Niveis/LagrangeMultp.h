#ifndef LAGRANGEMULTP_INCLUDED
#define LAGRANGEMULTP_INCLUDED  

#include "SolutionDecoder.h"

class LagrangeMultp: public SolutionDecoder{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
	int initInstance(InputFunction *function);
        int decodifySolution(Solution *sol);
};

#endif
 
