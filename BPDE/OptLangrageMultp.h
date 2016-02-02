#ifndef PENALTYFUNCTION_INCLUDED
#define PENALTYFUNCTION_INCLUDED  

#include "InputFunction.h"
#include "Solution.h"

double validateSolution(Solution *sol,InputFunction *function);
int compareSolutions(Solution *sol1, Solution *sol2,InputFunction *function);
int updateAPM(Solution **population, int sizePop, InputFunction *function);
int compareAPM(Solution *sol1, Solution *sol2,InputFunction *function);
#endif
 
