#include "InputFunction.h"
#include "Solution.h"
#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>

using namespace std;
int main(){
      InputFunction *function=new InputFunction("func1");
      srand(111111);
   
      DifferentialEvolution::initPopulation(function, 100);
      for(int i=0;i<5;i++){      
	DifferentialEvolution::mutatePopulation(0.3);
	DifferentialEvolution::recombinePopulation(0.8);
	DifferentialEvolution::selectPopulation(function);
      }
	
      if (DifferentialEvolution::best!=NULL){ 
	cout<<DifferentialEvolution::best->calcScore(function)<<"\n";
	for(int i=0;i<DifferentialEvolution::best->sizeVec;i++){
		cout<<DifferentialEvolution::best->vectorCharacters[i]<<"\t";
	}
      }else{
          cout<<"no feasible solution\n";
      }

      DifferentialEvolution::clearPopulation();
    

      delete function;
      
      
      
      return 1;
} 
