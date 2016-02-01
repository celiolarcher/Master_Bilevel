#include "InputFunction.h"
#include "Solution.h"
#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>


using namespace std;
int main(){
  InputFunction *function=new InputFunction("func1");
  srand(111111);
  
  /*
  double x[2]={19,14};
  double y[1]={14};
  
  cout<<function->getUPLevelFunction(x,y)<<"\n";

  cout<<function->getLWLevelFunction(x,y)<<"\n\n";

  double list[function->getNEQConstraintNumberLW()];
  function->constraintsValueNEQLW(x,y,list);

  for(int i=0;i<function->getNEQConstraintNumberLW();i++) cout<<list[i]<<"\n";
  
  Solution *sol=new Solution(function);
  
  sol->initValue(x,function);
  
  cout<<sol->calcScore(function)<<"\n";
  
  cout<<function->getLWLevelFunction(sol->vet,sol->vet + 1)<<"\n";
  */
 /* 
  Solution *best=NULL;
  for(int i=0;i<10000000000;i++){
        Solution *sol=new Solution(function);
        sol->initRandom(function);
        
        int flag=1;
        
    //    cout<<sol->vet[0]<<"\t"<<sol->vet[1]<<"\n";
        
        if(function->getNEQConstraintNumberLW()>0){
	double list[function->getNEQConstraintNumberLW()];
	function->constraintsValueNEQLW(sol->vet,sol->vet+function->getDimensionUP(),list);
	for(int i=0;i<function->getNEQConstraintNumberLW();i++){
	      if(list[i]>0){
	          flag=0;
	          break;
	      }
	}
        }
        
        if(function->getEQConstraintNumberLW()>0){
	double list[function->getEQConstraintNumberLW()];
	function->constraintsValueEQLW(sol->vet,sol->vet+function->getDimensionUP(),list);
	for(int i=0;i<function->getEQConstraintNumberLW();i++){
	      if(list[i]>0){
	          flag=0;
	          break;
	      }
	}
        }
        
        
        if(function->getNEQConstraintNumberUP()>0){
	double list[function->getNEQConstraintNumberUP()];
	function->constraintsValueNEQUP(sol->vet,sol->vet+function->getDimensionUP(),list);
	for(int i=0;i<function->getNEQConstraintNumberUP();i++){
	      if(list[i]>0){
	          flag=0;
	          break;
	      }
	}
        }
        
        
        if(function->getEQConstraintNumberUP()>0){
	double list[function->getEQConstraintNumberUP()];
	function->constraintsValueEQUP(sol->vet,sol->vet+function->getDimensionUP(),list);
	for(int i=0;i<function->getEQConstraintNumberUP();i++){
	      if(list[i]>0){
	          flag=0;
	          break;
	      }
	}
        }

        if(function->getKKTConstraintNumber()>0){
	double list[function->getKKTConstraintNumber()];
	function->constraintsValueKKT(sol->vet,sol->vet+function->getDimensionUP(),sol->vet+function->getDimensionUP()+function->getDimensionLW(),sol->vet+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
	for(int i=0;i<function->getKKTConstraintNumber();i++){
	      if(list[i]>-10e-2 && list[i]<10e-2){
	          flag=0;
	          break;
	      }
	}
        }

        if(function->getNEQConstraintNumberLW()>0){
	double list[function->getNEQConstraintNumberLW()];
	function->constraintsSlackness(sol->vet,sol->vet+function->getDimensionUP(),sol->vet+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
	
	for(int i=0;i<function->getNEQConstraintNumberLW();i++){
	      if(list[i]>-10e-1 && list[i]<10e-1){
	          flag=0;
	          break;
	      }
	}
        }
        
        if(flag) cout<<sol->calcScore(function)<<"\n";
        
        delete sol;
  }
  */
  
      DifferentialEvolution::initPopulation(function, 50);
      for(int i=0;i<5;i++){      
	DifferentialEvolution::mutatePopulation(0.5);
	DifferentialEvolution::recombinePopulation(0.6);
	DifferentialEvolution::selectPopulation(function);
      }
      DifferentialEvolution::clearPopulation();
    
      delete function;
      
      
      
      return 1;
} 
