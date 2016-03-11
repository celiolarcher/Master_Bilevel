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
	

/*
vectorCharacters[0]=1.95;
vectorCharacters[1]=8.05;
vectorCharacters[2]=0.00;
vectorCharacters[3]=0.97;
vectorCharacters[4]=0.97;
vectorCharacters[5]=1.31;
vectorCharacters[6]=6.74;
vectorCharacters[7]=0;
vectorCharacters[8]=0;
*/

//(0,2.44,10,0,10,8.74,5.25,10,0,10,3.73,10,10,10,0,0)

/*
vectorCharacters[0]=19;
vectorCharacters[1]=14;
vectorCharacters[2]=0;
vectorCharacters[3]=1;
vectorCharacters[4]=2;
vectorCharacters[5]=0;
*/
          return 1;
      }

      
      Solution *Solution::clone(){
           Solution *clone=new Solution(this->sizeVec,this->countConstraint);
	  for(int i=0;i<sizeVec;i++){
	      clone->vectorCharacters[i]=this->vectorCharacters[i];
	  }
              
	  for(int i=0;i<countConstraint;i++){
                     clone->constraintValues[i]=this->constraintValues[i];
               }

              clone->score=this->score;
              clone->upLevelFunction=this->upLevelFunction;
              clone->feasible=this->feasible;
              return clone;
      }
      
      std::ostream& operator<<(std::ostream &out, Solution &sol){

	  out<<"\n Solution Score:"<<sol.score<<"   Level UP:"<<sol.upLevelFunction<<"\n";
          
          for(int i=0;i<sol.sizeVec;i++) out<<sol.vectorCharacters[i]<<"\t";
          
          out<<"\n Constraints \n";          

          for(int i=0;i<sol.countConstraint;i++) out<<sol.constraintValues[i]<<"\t";
          
          out<<"\n";

          return out;
      }
