#ifndef DELOWERLEVEL_INCLUDED
#define DELOWERLEVEL_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"
#include "SolutionDecoder.h"
#include "DifferentialEvolution.h"

class DELowerLevel: public SolutionDecoder{
    private: DELowerLevel *UP;
	 DELowerLevel *LW;
	 DifferentialEvolution *DELW;
  
    public:
	int DEmiss;
	long long int DEsolution;
	int DECall;
	int initInstance(InputFunction *function);
	int decodifySolution(Solution *sol);
};

#endif
 
