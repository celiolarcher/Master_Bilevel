#ifndef LEMKELW_INCLUDED
#define LEMKELW_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"

class LemkeLW: public SolutionDecoder{
    public:
	int initInstance(InputFunction *function);
	int decodifySolution(Solution *sol);
	int getYLemke(double x[], double y[], double *opt);
};

#endif
 
