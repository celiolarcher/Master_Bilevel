#include "InputFunction.h"
#include "Solution.h"
#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>

using namespace std;
int main(){
      srand(111111);
      InputFunction *function=new InputFunction("func1");
      SolutionDecoder *decoder=new LagrangeMultpAPMSimplex();
      decoder->initInstance(function);
   
      DifferentialEvolution::initPopulation(decoder,50);
      for(int i=0;i<50000;i++){      
	DifferentialEvolution::mutatePopulation(0.6);
//	DifferentialEvolution::recombinePopulation(0.5);
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
