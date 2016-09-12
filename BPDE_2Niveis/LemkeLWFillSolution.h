#ifndef LEMKELWFILLSOLUTION_INCLUDED
#define LEMKELWFILLSOLUTION_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LemkeLWFillSolution: public SolutionDecoder{
    public:
	int initInstance(InputFunction *function);
	int decodifySolution(Solution *sol);
	int getYLemke(double x[], double y[], double *opt);
};

#endif
 
