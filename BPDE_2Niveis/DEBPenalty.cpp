#include "DEBPenalty.h"

#include <cmath>
#include <iostream>
using namespace std;

   
    int DEBPenalty::updatePenalty(Solution **population, Solution **newPopulation, int sizePop, int sizeNewPop, SolutionDecoder *decoder){
             
          return 1;
     }

     int DEBPenalty::compareSolutions(Solution *sol1, Solution *sol2, SolutionDecoder *decoder){   //1 se sol1 melhor que sol2, 0 caso contrário
          if(!sol1->completeSolution && sol2->completeSolution) return 0;
          if(sol1->completeSolution && !sol2->completeSolution) return 1;

          if(sol1->feasible && (!sol2->feasible || sol1->upLevelFunction < sol2->upLevelFunction)) return 1;
          else{
		  sol1->score=0;
		  sol2->score=0;

		  if(sol2->feasible) return 0;

		  for(int j=0;j<sol1->countConstraint;j++){
		          if(j < decoder->constraintsNEQNumber){
				  if(sol1->constraintValues[j]>0)
				      sol1->score+=sol1->constraintValues[j];
			  }else{
			      sol1->score+=fabs(sol1->constraintValues[j]);
		          }
		  }

		  for(int j=0;j<sol2->countConstraint;j++){
	          	  if(j < decoder->constraintsNEQNumber){
				  if(sol2->constraintValues[j]>0)
				      sol2->score+=sol2->constraintValues[j];
			  }else{
			      sol2->score+=fabs(sol2->constraintValues[j]);
		          }
		  }

		  if(!sol2->feasible && sol1->score < sol2->score) return 1;          
	  }
      
          return 0;
     }
