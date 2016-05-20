#include "DifferentialEvolution.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "DefaultParameters.h"

#define PENALTY_VALUE 1e4

void parseOptions(double *options, int dim[], int argc,char *argv[],char **file,long int *seed){ //Traduz os parâmetros recebidos pelo programa
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
        }else if(!strcmp(argv[i],"--mutImprov1")){
            options[3]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutImprov2")){
            options[4]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutFind1")){
            options[5]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutFind2")){
            options[6]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--crossRateImprov")){
            options[7]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutOptSearch")){
            options[8]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutBest1")){
            options[9]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutBest2")){
            options[10]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutBestSize")){
            options[11]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutOptBest")){
            options[12]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--mutOptImprov")){
            options[13]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--crossOptImprov")){
            options[14]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--crossOptSearch")){
            options[15]=atoi(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--crossRateSearch")){
            options[16]=atof(argv[i+1]);
            i+=2;
        }else if(!strcmp(argv[i],"--dim")){
            dim[0]=atoi(argv[i+1]);
	dim[1]=atoi(argv[i+2]);
	dim[2]=atoi(argv[i+3]);
	dim[3]=atoi(argv[i+4]);
            i+=5;
        }
        
        else i++;
    }
}


double fRand2(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}



extern double TOL_EQ_CONST;
extern double TOL_NEQ_CONST;
bool InfeasibleAvaliation;

using namespace std;
int main(int argc, char *argv[]){
      if(argc<2){
          srand(8990192031);
          InputFunction *function=new InputFunction("funcSMD1",3,2,2,0);
          //InputFunction *function=new InputFunction("funcJ2");
          SolutionDecoder *decoder=new LagrangeMultpSimplex();
          PenaltySolution *penalty=new APMDEBPenalty();
          decoder->initInstance(function);
      
          int popSize=80;
          DifferentialEvolution::initPopulation(decoder,penalty,popSize);
         // DifferentialEvolution::initPopulationNelderMeadMethod(decoder,penalty,popSize);   //FALTA FINALIZAR
          
          //DifferentialEvolution::improveInitSet(popSize,25, 0.53,0.23);
          DifferentialEvolution::improveInitSetSimilarity(popSize,25, 0.53,0.23,2,0.6,2,0);

          for(int i=0;i<50000 && function->getUPLevelCalls()<6000;i++){      
	//    DifferentialEvolution::mutatePopulationBestBounded(0.8,0.8);
	//    DifferentialEvolution::mutatePopulation_Rand_1_Bounded(0.8,0,50);
	    double rd=fRand2(0,0.2);
	    //DifferentialEvolution::mutatePopulationTargetBestBounded(0+rand,rand+0.6);



//	    if(function->getUPLevelCalls()<20)
//	      DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(0.8,0.8,0,popSize);
//	   else{
	    int ch=rand() % 2;
	      	//DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(0.9,0.8,0,20);
	    //if(ch)
	      DifferentialEvolution::mutatePopulation_RandToBest_1_Bounded(rd+0.1,1.0+rd,0,25);
	    //else
	      //DifferentialEvolution::mutatePopulation_Best_2_Bounded(0.7,0.7,0,25);
	
	//	DifferentialEvolution::mutatePopulation_BestToRand_1_Bounded(0.8,0.4,15,20);
//	    }
/*
	    if(function->getUPLevelCalls()>0){
//		cout<<"aqui";
		popSize=20;
		DifferentialEvolution::sizePopulation=20;
	    }
*/
	//    DifferentialEvolution::recombinePopulation(0.9);	   
//	    	DifferentialEvolution::recombinePopulationExp();
	    	DifferentialEvolution::selectPopulation();
          }
	    
          if (DifferentialEvolution::best!=NULL){ 
	  	cout<<*DifferentialEvolution::best<<"UP Level Calls Until Best:"<<DifferentialEvolution::UPLevelCallsBest<<"\nUP Level Calls Process: "<<function->getUPLevelCalls()<<"\n";
		cout<<"UP Level \t"<<function->getUPLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";
		cout<<"Lower Level \t"<<function->getLWLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";
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

          long int seed=911123193131;
          char *file;
          int dim[4]={-1,-1,-1,-1};
          
          
          parseOptions(options,dim,argc,argv,&file,&seed);
          for(int i=0;i<COUNTPARAMETERS;i++){
	        if(options[i]==-1)options[i]=defaultParameters[i];
          }

          srand(seed);
          
          InputFunction *function;
          if(dim[0]<0)
	function=new InputFunction(file);
          else
	function=new InputFunction(file,dim[0],dim[1],dim[2],dim[3]);
          

          SolutionDecoder *decoder=new LemkeLW();
          PenaltySolution *penalty=new DEBPenalty();
          decoder->initInstance(function);

          TOL_EQ_CONST=1e-4;
          TOL_NEQ_CONST=1e-4;
          
          InfeasibleAvaliation=0;
          
          //DifferentialEvolution::initPopulation(decoder,penalty,options[0]);  
          
          DifferentialEvolution::initPopulation(decoder,penalty,options[1]);  //Versão sem aprimoramento inicial
          
          //DifferentialEvolution::initPopulationNelderMeadMethod(decoder,penalty,options[0]);
          
          
       //   DifferentialEvolution::improveInitSetSimilarity(options[0],options[1],options[5],options[6],options[8],options[16],options[15],options[11]);          
        //  DifferentialEvolution::improveInitSetDispersion(options[0],options[1],options[5],options[6],options[8],options[16],options[15],options[11]);

         // TOL_EQ_CONST=1e-4;
//          // TOL_NEQ_CONST=1e-4;
          
         // Solution *bestTurn=DifferentialEvolution::best->clone();
                    
//          DifferentialEvolution::decodifyPopulation();
          
          for(int i=0;i<100000 && function->getUPLevelCalls()<options[2] && function->getLWLevelSimplexCalls()<1e7;i++){      
	  /*
	          if(function->getUPLevelCalls()>1500){
		  //TOL_EQ_CONST=1e-4;
		  TOL_NEQ_CONST=1e-4;
		  DifferentialEvolution::decodifyPopulation();
	          }
	
	/*
	          if(function->getUPLevelCalls()>2500){
		  TOL_EQ_CONST=1e-4;
		  TOL_NEQ_CONST=1e-4;
		  DifferentialEvolution::decodifyPopulation();
	          }*/
	          
	
		double rd=fRand2(0,options[4]);	

	
		if(options[13]==1)
		  DifferentialEvolution::mutatePopulation_TargetToBest_1(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==2)
		  DifferentialEvolution::mutatePopulation_TargetToRand_1(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==3)
		  DifferentialEvolution::mutatePopulation_Target_1(options[3]+rd,0,options[1]);
		else if(options[13]==4)
		  DifferentialEvolution::mutatePopulation_RandToBest_1(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==5)
		  DifferentialEvolution::mutatePopulation_Rand_1(options[3]+rd,0,options[1]);
		else if(options[13]==6)
		  DifferentialEvolution::mutatePopulation_Rand_2(options[3]+rd,0,options[1]);
		else if(options[13]==7)  
		  DifferentialEvolution::mutatePopulation_Target_2(options[3]+rd,options[3]+rd,0,options[1]);
		
		//versão bounded 
		else if(options[13]==8)
		  DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(options[3]+rd,options[3]+rd,0,options[1]); 
		else if(options[13]==9)
		  DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==10)
		  DifferentialEvolution::mutatePopulation_Target_1_Bounded(options[3]+rd,0,options[1]);
		else if(options[13]==11)
		  DifferentialEvolution::mutatePopulation_RandToBest_1_Bounded(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==12)
		  DifferentialEvolution::mutatePopulation_Rand_1_Bounded(options[3]+rd,0,options[1]);
		else if(options[13]==13)
		  DifferentialEvolution::mutatePopulation_Rand_2_Bounded(options[3]+rd,0,options[1]);
		else if(options[13]==14)  
		  DifferentialEvolution::mutatePopulation_Target_2_Bounded(options[3]+rd,options[3]+rd,0,options[1]);
		
	/*
		if(options[13]==1)
		  DifferentialEvolution::mutatePopulation_RandToBest_1(options[3]+rd,options[3]+rd,0,options[1]);
		else if(options[13]==2)
		  DifferentialEvolution::mutatePopulation_Rand_1(options[3]+rd,0,options[1]);
		else if(options[13]==3)
		  DifferentialEvolution::mutatePopulation_Rand_2(options[3]+rd,0,options[1]);
		else if(options[13]==4)
		  DifferentialEvolution::mutatePopulation_RandToBest_1_Bounded(options[3]+rd,options[3],0,options[1]);
		else if(options[13]==5)
		  DifferentialEvolution::mutatePopulation_Rand_1_Bounded(options[3]+rd,0,options[1]);
		else if(options[13]==6)
		  DifferentialEvolution::mutatePopulation_Rand_2_Bounded(options[3]+rd,0,options[1]);
	*/
	
	
		if(options[12]==1)
		  DifferentialEvolution::mutatePopulation_BestToRand_1(options[9],options[10],options[1],options[1]+options[11]);
		if(options[12]==2)
		  DifferentialEvolution::mutatePopulation_Best_1(options[9],options[1],options[1]+options[11]);
		if(options[12]==3)
		  DifferentialEvolution::mutatePopulation_Best_2(options[9],options[10],options[1],options[1]+options[11]);

		if(options[14]==1)
		      DifferentialEvolution::recombinePopulationSwap(0,options[1]);
		else if(options[14]==2)
		      DifferentialEvolution::recombinePopulation(options[7],0,options[1]);
		else if(options[14]==4)
		      DifferentialEvolution::recombinePopulationExp(options[7],0,options[1]);
		
		DifferentialEvolution::selectPopulation();
		/*
		if(bestTurn->diffMaxSolution(DifferentialEvolution::best)>1e-4){
		    cout<<"it "<<i<<"\n";
		      delete bestTurn;
		      bestTurn=DifferentialEvolution::best->clone();
		      
		      for(int i=options[1];i<options[1]+options[11];i++){
		          cout<<"----------------------------------------\n"<<*DifferentialEvolution::Population[i];
		          delete DifferentialEvolution::Population[i];
		          DifferentialEvolution::Population[i]=DifferentialEvolution::best->clone();
		          cout<<*DifferentialEvolution::Population[i];
		      }
		      
		
		    cout<<"\n\n";
		  
		}*/
		
          }
	    

          if (DifferentialEvolution::best!=NULL && DifferentialEvolution::best->feasible){ 
		    cout<<DifferentialEvolution::best->upLevelFunction<<"\t"<<DifferentialEvolution::UPLevelCallsBest<<"\t"<<DifferentialEvolution::LWLevelSimplexCallsBest<<"\t"<<function->getLWLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP());
		    cout<<"\t(";
		    for(int j=0;j<DifferentialEvolution::best->sizeVec-1;j++)
		        cout<<DifferentialEvolution::best->vectorCharacters[j]<<",";
		    cout<<DifferentialEvolution::best->vectorCharacters[DifferentialEvolution::best->sizeVec-1]<<")";
		    
		    		//cout<<*DifferentialEvolution::best;
		//		cout<<"UP Level \t"<<function->getUPLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";
		//cout<<"Lower Level \t"<<function->getLWLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";

          }else{
		cout<<PENALTY_VALUE;
		//cout<<*DifferentialEvolution::best;
		//		cout<<"UP Level \t"<<function->getUPLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";
		//cout<<"Lower Level \t"<<function->getLWLevelFunction(DifferentialEvolution::best->vectorCharacters,DifferentialEvolution::best->vectorCharacters+function->getDimensionUP())<<"\n";

//		  cout<<"no feasible solution";
          }

          DifferentialEvolution::clearPopulation();
        
          delete function;
          delete decoder;
          delete penalty;
        
      }
      
      
      return 1;
} 
