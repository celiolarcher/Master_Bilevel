#include "InputFunction.h"
#include "Solution.h"
#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>

using namespace std;
int main(){
      srand(111111);
      InputFunction *function=new InputFunction("func1");
      SolutionDecoder *decoder=new LagrangeMultpAPMSmooth();
      decoder->initInstance(function);
   
      DifferentialEvolution::initPopulation(decoder,100);
      for(int i=0;i<1;i++){      
	DifferentialEvolution::mutatePopulation(0.5);
	DifferentialEvolution::recombinePopulation(0.6);
	DifferentialEvolution::selectPopulation();
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
