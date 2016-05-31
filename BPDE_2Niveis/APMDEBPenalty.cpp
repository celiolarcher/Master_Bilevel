#include "APMDEBPenalty.h"

#include <cmath>
#include <iostream>
using namespace std;

  
    int APMDEBPenalty::updatePenalty(Solution **population, Solution **newPopulation, int sizePop, int sizeNewPop, SolutionDecoder *decoder){
          levelUPMean=0;
          
          double violationMean[population[0]->countConstraint];
          if(!kWeight)
		kWeight=new double[population[0]->countConstraint];
          
          for(int j=0;j<population[0]->countConstraint;violationMean[j]=0, j++);
	  
          for(int i=0;i<sizePop;i++){
		  levelUPMean+=population[i]->upLevelFunction;
		  for(int j=0;j<population[0]->countConstraint;j++){
			if(j < decoder->constraintsNEQNumber){
				if(population[i]->constraintValues[j]>0){
				    violationMean[j]+=population[i]->constraintValues[j];
				} 
			}else{
				violationMean[j]+=fabs(population[i]->constraintValues[j]);
			}
		  }
          }
          
          if(newPopulation){
		  for(int i=0;i<sizeNewPop;i++){
		      levelUPMean+=newPopulation[i]->upLevelFunction;
		      for(int j=0;j<newPopulation[0]->countConstraint;j++){
			if(j < decoder->constraintsNEQNumber){
				    if(newPopulation[i]->constraintValues[j]>0){
					violationMean[j]+=newPopulation[i]->constraintValues[j];
				    } 
			}else{
				violationMean[j]+=fabs(newPopulation[i]->constraintValues[j]);
			}
		      }
		  }
          }
          
//          cout<<"\n--------------------------------------------------------------------------\n";
          double violationSumSquare=0;
          for(int j=0;j<population[0]->countConstraint;j++){
		if(newPopulation)
		    violationMean[j]/=(2*sizePop);
		else
		    violationMean[j]/=(sizePop);
	   
	        violationSumSquare+=violationMean[j]*violationMean[j];
  //        cout<<"\t"<<violationSumSquare;                    
          }
//cout<<"\n";
          
        if(newPopulation)
          levelUPMean/=(2*sizePop);
	else
	  levelUPMean/=(sizePop);
          
          for(int j=0;j<population[0]->countConstraint;j++){
//	     kWeight[j]=(violationMean[j]*fabs(levelUPMean))/violationSumSquare;
	     kWeight[j]=(violationMean[j])/violationSumSquare;
		//cout<<kWeight[j]<<"\t";
          }
          
          //cout<<endl;
          //cout<<levelUPMean<<"\t"<<violationSumSquare<<"\n";          
          return 1;
     }

     int APMDEBPenalty::compareSolutions(Solution *sol1, Solution *sol2, SolutionDecoder *decoder){   //1 se sol1 melhor que sol2, 0 caso contrário
	  if(!sol1->completeSolution) return 0;
	  if(!sol2->completeSolution) return 1;

          if(sol1->feasible) sol1->score=sol1->upLevelFunction;
          else{
	    sol1->score=0;
	    for(int j=0;j<sol1->countConstraint;j++){
	          if(j < decoder->constraintsNEQNumber){
		  if(sol1->constraintValues[j]>0)
		      sol1->score+=kWeight[j]*sol1->constraintValues[j];
	          }else{
		      sol1->score+=kWeight[j]*fabs(sol1->constraintValues[j]);
	          }
	    }

	    //if(sol1->upLevelFunction>levelUPMean)sol1->score+=sol1->upLevelFunction;
	    //else sol1->score+=levelUPMean;
	    sol1->score+=sol1->penaltyValue;
          }
          
      
          if(sol2->feasible) sol2->score=sol2->upLevelFunction;
          else{
	    sol2->score=0;
	    for(int j=0;j<sol2->countConstraint;j++){
	          if(j < decoder->constraintsNEQNumber){
		  if(sol2->constraintValues[j]>0)
		      sol2->score+=kWeight[j]*sol2->constraintValues[j];
	          }else{
		      sol2->score+=kWeight[j]*fabs(sol2->constraintValues[j]);
	          }
	    }

	    //if(sol2->upLevelFunction>levelUPMean)sol2->score+=sol2->upLevelFunction;
	    //else sol2->score+=levelUPMean;
	    sol2->score+=sol2->penaltyValue;
          }

	  if(sol1->feasible && !sol2->feasible) return 1;
	  if(!sol1->feasible && sol2->feasible) return 0;

	
          if(sol1->score < sol2->score) return 1;
          
          return 0;
     }
