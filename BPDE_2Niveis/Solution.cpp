#include "Solution.h"
#include <stdlib.h>
#include <cmath>

extern double TOL_EQ_CONST=1e-4;
extern double TOL_NEQ_CONST=1e-4;

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
          penaltyValue=0;
          return 1;
      }
      
      int Solution::initRandom(double bounds[]){
          for(int i=0;i<sizeVec;i++){
	  vectorCharacters[i]=fRand(bounds[2*i],bounds[2*i+1]);
          }
          
          
          //vectorCharacters[0]=0.5;
      //    ,,,,0,10.6102,34.2204,1.22037
          
           //  vectorCharacters[4]=0;
          
          
          /*
          vectorCharacters[0]=3;
          vectorCharacters[1]=3;
          vectorCharacters[2]=3;
      //     vectorCharacters[0]=3;
  /*
          vectorCharacters[0]=0.452808;
          vectorCharacters[1]=0.40636;
          vectorCharacters[2]=9.63507;
          vectorCharacters[3]=0.01392;
          vectorCharacters[4]=9.62655;
                 vectorCharacters[5]=8.92148;
	      vectorCharacters[6]=2.66284;
vectorCharacters[7]=9.92573;
	      vectorCharacters[8]=0.205161;
	      vectorCharacters[9]=9.84991;

  
        //(,,,,,,,,,,10,10,9.6916,9.99554,-2.11738e-11,7.93686)
/*        vectorCharacters[0]=0.452808;
          vectorCharacters[1]=0.40636;
          vectorCharacters[2]=9.63507;
          vectorCharacters[3]=0.01392;
          vectorCharacters[4]=9.62655;
                 vectorCharacters[5]=8.92148;
	      vectorCharacters[6]=2.66284;
vectorCharacters[7]=9.92573;
	      vectorCharacters[8]=0.205161;
	      vectorCharacters[9]=9.84991;

          
          /*vectorCharacters[0]=0.5;
        vectorCharacters[1]=0.8;
       
     /*     vectorCharacters[0]=0;
          vectorCharacters[1]=2.44;
          vectorCharacters[2]=10;
          vectorCharacters[3]=0;
          vectorCharacters[4]=10;
          vectorCharacters[5]=8.74;
          vectorCharacters[6]=5.25;
          vectorCharacters[7]=10;
          vectorCharacters[8]=0;
          vectorCharacters[9]=10;

/*          
vectorCharacters[0]=1.95;
vectorCharacters[1]=8.05;
vectorCharacters[2]=0;
vectorCharacters[3]=0.97;
vectorCharacters[4]=0.97;
vectorCharacters[5]=1.31;
vectorCharacters[6]=6.74;
vectorCharacters[7]=0;
vectorCharacters[8]=0;
*/
  //       vectorCharacters[0]=0.5;
//vectorCharacters[1]=0.8;

         
/*
         vectorCharacters[0]=0.5;
vectorCharacters[1]=0.8;
vectorCharacters[2]=0;
vectorCharacters[3]=0.2;
vectorCharacters[4]=0.8;
         */
         
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
          penaltyValue=0;

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
              clone->completeSolution=this->completeSolution;
              return clone;
      }
      
      std::ostream& operator<<(std::ostream &out, Solution &sol){

	  out<<"\n Solution Score:"<<sol.score<<"   Level UP:"<<sol.upLevelFunction<<"    Complete Solution:"<<sol.completeSolution<<"    Feasible Solution:"<<sol.feasible<<"\n";
          
          for(int i=0;i<sol.sizeVec;i++) out<<sol.vectorCharacters[i]<<"\t";
          
          out<<"\n Constraints \n";          

          for(int i=0;i<sol.countConstraint;i++) out<<sol.constraintValues[i]<<"\t";
          
          out<<"\n";

          return out;
      }

      
      double Solution::diffZeroSolution(Solution *sol){
          double diff=0;
          
          for(int i=0;i<sizeVec;i++){
	  diff+=fabs(this->vectorCharacters[i]-sol->vectorCharacters[i]);
          }
          
          return diff;
      }
      
      double Solution::diffSquareSolution(Solution *sol){
          double diff=0;
          
          for(int i=0;i<sizeVec;i++){
	  diff+=(this->vectorCharacters[i]-sol->vectorCharacters[i])*(this->vectorCharacters[i]-sol->vectorCharacters[i]);
          }
          
          return diff;
      }
      
      
            
      double Solution::diffMaxSolution(Solution *sol){
          double diff=0;
          
          for(int i=0;i<sizeVec;i++){
	
	  double aux=fabs(this->vectorCharacters[i]-sol->vectorCharacters[i]);
	  
	  if(diff<aux)  diff=aux;
          }
          
          return diff;
      }
