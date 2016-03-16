#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>

double fRand2(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

using namespace std;
int main(int argc, char *argv[]){
      if(argc<2){
          srand(31231231);
          InputFunction *function=new InputFunction("funcB8");
          SolutionDecoder *decoder=new LagrangeMultpSimplex();
          PenaltySolution *penalty=new APMDEBPenalty();
          decoder->initInstance(function);
      
          DifferentialEvolution::initPopulation(decoder,penalty,50);
          for(int i=0;i<500000 && function->getUPLevelCalls()<5000;i++){      
	//    DifferentialEvolution::mutatePopulationBestBounded(0.8,0.8);
	    //DifferentialEvolution::mutatePopulationTargetBounded(0.8,0.8);
	    double rand=fRand2(0,0.2);
	    //DifferentialEvolution::mutatePopulationTargetBestBounded(0+rand,rand+0.6);
	    if(function->getUPLevelCalls()<2500)
	      DifferentialEvolution::mutatePopulationTargetBounded(0.8,0.8);
	    else
	      DifferentialEvolution::mutatePopulationBestBounded(0.8,0.8);
	    //DifferentialEvolution::mutatePopulationBounded(0.8);
	   
	    DifferentialEvolution::recombinePopulation(0.6);
	    DifferentialEvolution::selectPopulation();
          }
	    
          if (DifferentialEvolution::best!=NULL){ 
	  cout<<*DifferentialEvolution::best<<"UP Level Calls:"<<function->getUPLevelCalls()<<"\n";
          }else{
	  cout<<"no feasible solution\n";
          }

          DifferentialEvolution::clearPopulation();
        
          delete function;
          delete decoder;
          delete penalty;
      }else{
          srand(atoi(argv[1]));
          InputFunction *function=new InputFunction(argv[2]);
          SolutionDecoder *decoder=new LagrangeMultpSimplex();
          PenaltySolution *penalty=new APMDEBPenalty();
          decoder->initInstance(function);
      
          DifferentialEvolution::initPopulation(decoder,penalty,50);
          for(int i=0;i<500000 && function->getUPLevelCalls()<5000;i++){      
	    //DifferentialEvolution::mutatePopulationTargetBestBounded(0.3,0.8);
	    //DifferentialEvolution::mutatePopulationBestBounded(0.3,0.8);
	    //DifferentialEvolution::mutatePopulationBounded(0.7);
	    //DifferentialEvolution::mutatePopulationTargetBounded(0.8,0.8);
	     if(function->getUPLevelCalls()<2500)
	      DifferentialEvolution::mutatePopulationTargetBounded(0.8,0.8);
	    else
	      DifferentialEvolution::mutatePopulationBestBounded(0.8,0.8);
	    DifferentialEvolution::recombinePopulation(0.6);
	    DifferentialEvolution::selectPopulation();
          }
	    
          if (DifferentialEvolution::best!=NULL && DifferentialEvolution::best->feasible){ 
		    cout<<DifferentialEvolution::best->upLevelFunction;
          }else{
	  cout<<"no feasible solution";
          }

          DifferentialEvolution::clearPopulation();
        
          delete function;
          delete decoder;
          delete penalty;
        
      }
      
      
      return 1;
} 
