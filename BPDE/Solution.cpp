#include "Solution.h"
#include <stdlib.h>

      double fRand(double fMin, double fMax){
          double f = (double)rand() / RAND_MAX;
          return fMin + f * (fMax - fMin);
      }


      double Solution::calcScore(InputFunction *function){
          return function->getUPLevelFunction(vet,vet+function->getDimensionUP());
      }

      int Solution::initValue(double initVet[], InputFunction *function){
          for(int i=0;i<function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();i++){
	  	vet[i]=initVet[i];
          }
          return 1;
      }
      
      int Solution::initRandom(InputFunction *function){
        
          for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	 	 vet[i]=fRand(10,20);
          }

          
          for(int i=function->getDimensionUP()+function->getDimensionLW();i<function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();i++){
	  vet[i]=fRand(0,2);
          }
          
          return 1;
      }
      
      
