#ifndef LAGRANGEMULTPSIMPLEX_INCLUDED
#define LAGRANGEMULTPSIMPLEX_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LagrangeMultpSimplex: public SolutionDecoder{
    public:
    	double levelUPMean;
    	double *kWeight=NULL;
    public:
	int initInstance(InputFunction *function);
        int decodifySolution(Solution *sol);
	int simplexLagrangeMultipliers(double x[], double y[], double h[], double *opt);
};

#endif
 
