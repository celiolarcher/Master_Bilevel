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
   
      DifferentialEvolution::initPopulation(decoder,30);
      for(int i=0;i<500;i++){      
	DifferentialEvolution::mutatePopulation(0.8);
	DifferentialEvolution::recombinePopulation(0.1);
	DifferentialEvolution::selectPopulation();
      }
	
      if (DifferentialEvolution::best!=NULL){ 
		cout<<*DifferentialEvolution::best;
      }else{
          cout<<"no feasible solution\n";
      }

      DifferentialEvolution::clearPopulation();
    

      delete function;
      
      
      
      return 1;
} 
