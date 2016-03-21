#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "DefaultParameters.h"

#define PENALTY_VALUE 10e3

void parseOptions(double *options, int argc,char *argv[],char **file,long int *seed){ //Traduz os parâmetros recebidos pelo programa
    int i=0;
    while(i<argc){
         if(!strcmp(argv[i],"-i")){
            (*file)=argv[i+1];
	    char *find=strchr((*file),'/');
	    while (find!=NULL){
		(*file)=find+1;
	    	find=strchr(find+1,'/');
	    }

            i+=2;
        }
	else if(!strcmp(argv[i],"--seed")){
            (*seed)=atoi(argv[i+1]);
            i+=2;
	}
        else if(!strcmp(argv[i],"--populationSearch")){
            options[0]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--populationImprov")){
            options[1]=atoi(argv[i+1]);
            i+=2;
        }
        else if(!strcmp(argv[i],"--maxUPCalls")){
            options[2]=atoi(argv[i+1]);
            i+=2;
        }
	else if(!strcmp(argv[i],"--mutTargetRate")){
            options[3]=atof(argv[i+1]);
            i+=2;
        }
	else if(!strcmp(argv[i],"--mutBestRate")){
            options[4]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutFind1")){
            options[5]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutFind2")){
            options[6]=atof(argv[i+1]);
            i+=2;
        }
        else if(!strcmp(argv[i],"--crossRate")){
            options[7]=atof(argv[i+1]);
            i+=2;
        }
        else i++;
    }
}


double fRand2(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

using namespace std;
int main(int argc, char *argv[]){
      if(argc<2){
          srand(9429492);
          InputFunction *function=new InputFunction("funcB4");
          SolutionDecoder *decoder=new LagrangeMultpSimplex();
          PenaltySolution *penalty=new APMDEBPenalty();
          decoder->initInstance(function);
      
	  int popSize=100;
          DifferentialEvolution::initPopulation(decoder,penalty,popSize);
   	  DifferentialEvolution::improveInitSet(100,20, 0.7,0.7);

          for(int i=0;i<500000 && function->getUPLevelCalls()<20000;i++){      
	//    DifferentialEvolution::mutatePopulationBestBounded(0.8,0.8);
	//    DifferentialEvolution::mutatePopulation_Rand_1_Bounded(0.8,0,50);
	  //  double rand=fRand2(0,0.2);
	    //DifferentialEvolution::mutatePopulationTargetBestBounded(0+rand,rand+0.6);



//	    if(function->getUPLevelCalls()<20)
//	      DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(0.8,0.8,0,popSize);
//	   else{
	      	DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(0.9,0.8,0,15);
		DifferentialEvolution::mutatePopulation_BestToRand_1_Bounded(0.8,0.4,15,20);
//	    }
/*
	    if(function->getUPLevelCalls()>0){
//		cout<<"aqui";
		popSize=20;
		DifferentialEvolution::sizePopulation=20;
	    }
*/
//	    DifferentialEvolution::recombinePopulation(0.6);	   
	    	DifferentialEvolution::recombinePopulationExp();
	    	DifferentialEvolution::selectPopulation();
          }
	    
          if (DifferentialEvolution::best!=NULL){ 
	  	cout<<*DifferentialEvolution::best<<"UP Level Calls Until Best:"<<DifferentialEvolution::UPLevelCallsBest<<"\nUP Level Calls Process: "<<function->getUPLevelCalls()<<"\n";
          }else{
	  	cout<<"no feasible solution\n";
          }

          DifferentialEvolution::clearPopulation();
        
          delete function;
          delete decoder;
          delete penalty;




      }else{




	  double *options=new double[COUNTPARAMETERS];
	  for(int i=0;i<COUNTPARAMETERS;i++)options[i]=-1;

	  long int seed=3182391883;
	  char *file;

       	  parseOptions(options,argc,argv,&file,&seed);
	  for(int i=0;i<COUNTPARAMETERS;i++){
		if(options[i]==-1)options[i]=defaultParameters[i];
	  }

          srand(seed);
          InputFunction *function=new InputFunction(file);


          SolutionDecoder *decoder=new LagrangeMultpSimplex();
          PenaltySolution *penalty=new APMDEBPenalty();
          decoder->initInstance(function);



      
          DifferentialEvolution::initPopulation(decoder,penalty,options[0]);
   	  DifferentialEvolution::improveInitSet(options[0],options[1],options[5],options[6]);

          for(int i=0;i<100000 && function->getUPLevelCalls()<options[2];i++){      

	      	DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(options[4],options[3],0,2*options[1]/3);
		DifferentialEvolution::mutatePopulation_BestToRand_1_Bounded(options[4],options[3],2*options[1]/3,options[1]);

	    	DifferentialEvolution::recombinePopulationExp();
	    	DifferentialEvolution::selectPopulation();
          }
	    

          if (DifferentialEvolution::best!=NULL && DifferentialEvolution::best->feasible){ 
		    cout<<DifferentialEvolution::best->upLevelFunction<<"\t"<<DifferentialEvolution::UPLevelCallsBest;
          }else{
		cout<<PENALTY_VALUE;
//		  cout<<"no feasible solution";
          }

          DifferentialEvolution::clearPopulation();
        
          delete function;
          delete decoder;
          delete penalty;
        
      }
      
      
      return 1;
} 
