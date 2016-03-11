#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>

using namespace std;
int main(){
      srand(111111);
      InputFunction *function=new InputFunction("funcB5");
      SolutionDecoder *decoder=new LagrangeMultpSimplex();
      PenaltySolution *penalty=new APMDEBPenalty();
      decoder->initInstance(function);
   
      DifferentialEvolution::initPopulation(decoder,penalty,50);
      for(int i=0;i<50000;i++){      
	DifferentialEvolution::mutatePopulationBounded(0.8);
	DifferentialEvolution::recombinePopulation(0.6);
	DifferentialEvolution::selectPopulation();
      }
	
      if (DifferentialEvolution::best!=NULL){ 
		cout<<*DifferentialEvolution::best;
      }else{
          cout<<"no feasible solution\n";
      }

      DifferentialEvolution::clearPopulation();
    
      delete function;
      delete decoder;
      delete penalty;
      
      
      
      return 1;
} 
