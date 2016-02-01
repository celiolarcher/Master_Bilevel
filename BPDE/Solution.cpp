#include "Solution.h"


#include <stdlib.h>
      double fRand(double fMin, double fMax){
          double f = (double)rand() / RAND_MAX;
          return fMin + f * (fMax - fMin);
      }


      double Solution::calcScore(InputFunction *function){
          return function->getUPLevelFunction(vectorCharacters,vectorCharacters+function->getDimensionUP());
      }

      int Solution::initValue(double initVec[]){
          for(int i=0;i<sizeVec;i++){
	  vectorCharacters[i]=initVec[i];
          }
          return 1;
      }
      
      int Solution::initRandom(double bounds[]){
        
          for(int i=0;i<sizeVec;i++){
	  vectorCharacters[i]=fRand(bounds[2*i],bounds[2*i+1]);
          }
          
          return 1;
      }
      
      
      std::ostream& operator<<(std::ostream &out, Solution &sol){
          
          for(int i=0;i<sol.sizeVec;i++) out<<sol.vectorCharacters[i]<<"\t";
          
          out<<"\n";
          
          return out;
      }
