#include "DifferentialEvolution.h" 
#include "stdlib.h"
#define TRUNCCASE 10e5

#include <iostream>
#include <cmath>
using namespace std;
      Solution **DifferentialEvolution::Population;
      Solution **DifferentialEvolution::nextPopulation;
      Solution *DifferentialEvolution::best=NULL;
      int DifferentialEvolution::UPLevelCallsBest;
      int DifferentialEvolution::sizePopulation;
      SolutionDecoder *DifferentialEvolution::decoder;
      PenaltySolution *DifferentialEvolution::penalty;


      int DifferentialEvolution::initPopulation(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop){
	Population=new Solution*[sizePop];
	nextPopulation=new Solution*[sizePop];
	sizePopulation=sizePop;
	DifferentialEvolution::decoder=decoder;
	DifferentialEvolution::penalty=penalty;
	
	for(int i=0;i<sizePopulation; Population[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), nextPopulation[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), Population[i++]->initRandom(decoder->boundAttributes)); 
        
	
	for(int i=0;i<sizePopulation;i++){
	    decoder->decodifySolution(Population[i]);
	}
	
	penalty->updatePenalty(Population,NULL,sizePopulation,0,decoder);
	for(int i=0;i<sizePopulation;i++){
		//cout<<*Population[i];
	    if((best==NULL || penalty->compareSolutions(Population[i],best,decoder))){
	        if(!best) delete best;
	        best=Population[i]->clone();
		UPLevelCallsBest=decoder->function->getUPLevelCalls();
	    }
	}
	
	return 1;
      }
      
      


     int DifferentialEvolution::improveInitSet(int sizePopSearch, int sizePopNextStep, double find1, double find2){
	
		Solution **populationNew=new Solution *[sizePopNextStep];
		int feasibleSolutions=0;		

		for(int i=0;i<100000 && feasibleSolutions<sizePopNextStep/2; i++){
			DifferentialEvolution::mutatePopulation_Target_2_Bounded(find1,find2,0,sizePopSearch-feasibleSolutions);
			DifferentialEvolution::recombinePopulationExp();
			DifferentialEvolution::selectPopulation();

			for(int j=0;j<sizePopSearch-feasibleSolutions && feasibleSolutions<sizePopNextStep;j++){
				if(Population[j]->feasible){
					//cout<<j<<"\t"<<feasibleSolutions<<"\t"<<sizePopSearch-feasibleSolutions-1<<"\n";
					populationNew[feasibleSolutions]=Population[j];
					Population[j]=Population[sizePopSearch-feasibleSolutions-1];
					sizePopulation--;
					feasibleSolutions++;
				}
			}
		}
	//	cout<<feasibleSolutions;

		for(int j=feasibleSolutions;j<sizePopNextStep;j++) populationNew[j]=Population[j-feasibleSolutions];
	
		for(int j=sizePopNextStep-feasibleSolutions;j<sizePopSearch-sizePopNextStep;j++) delete Population[j];

		Solution **nextPopulationNew=new Solution*[sizePopNextStep];


		for(int j=0;j<sizePopNextStep;j++)nextPopulationNew[j]=nextPopulation[j];


		for(int j=sizePopSearch;j<sizePopSearch;j++) delete nextPopulation[j];

		delete Population;
		delete nextPopulation;
		
		Population=populationNew;
		nextPopulation=nextPopulationNew;
		sizePopulation=sizePopNextStep;
	
		return 1;
      }
      









      int DifferentialEvolution::mutatePopulation_Rand_1(double F, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + F*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		//if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) nextPopulation[i]->vectorCharacters[j]=decoder->boundAttributes[2*j];//Population[r1]->vectorCharacters[j];
		//if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) nextPopulation[i]->vectorCharacters[j]=decoder->boundAttributes[2*j+1];//Population[r1]->vectorCharacters[j];

	    }

	}
	
	return 1;
      }
      
      int DifferentialEvolution::mutatePopulation_Rand_1_Bounded(double F,int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=F;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	            nextPopulation[i]->vectorCharacters[j]=(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (Population[r1]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-Population[r1]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (Population[r1]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-Population[r1]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	//	if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) cout<<Population[r1]->vectorCharacters[j]<<"\t"<<minF<<"\t"<<(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j])<<"\t"<<nextPopulation[i]->vectorCharacters[j]<<"\n";
	//	if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) cout<<Population[r1]->vectorCharacters[j]<<"\t"<<minF<<"\t"<<(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j])<<"\t"<<nextPopulation[i]->vectorCharacters[j]<<"\n";

//		if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) cout<<nextPopulation[i]->vectorCharacters[j]-decoder->boundAttributes[2*j]<<"\n";
//		if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) cout<<nextPopulation[i]->vectorCharacters[j]-decoder->boundAttributes[2*j]<<"\n";

	    }
//if(minF<0)cout<<*nextPopulation[i];	

	}
	
	return 1;
      }
      

      int DifferentialEvolution::mutatePopulation_RandToBest_1_Bounded(double F1, double F2,int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
		nextPopulation[i]->vectorCharacters[j]=F1*(best->vectorCharacters[j] - Population[r1]->vectorCharacters[j])+F2*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (Population[r1]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-Population[r1]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
			//cout<<minF<<"\t";
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (Population[r1]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-Population[r1]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }

	}
	
	return 1;
      }
      
      int DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(double F1, double F2, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[i]->vectorCharacters[j])+F2*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
			//cout<<minF<<"\t";
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }


	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(double F1, double F2, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(best->vectorCharacters[j] - Population[i]->vectorCharacters[j])+F2*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
			//cout<<minF<<"\t";
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }


      int DifferentialEvolution::mutatePopulation_BestToRand_1_Bounded(double F1, double F2, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j]-best->vectorCharacters[j])+F2*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (best->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-best->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
			//cout<<minF<<"\t";
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (best->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-best->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=best->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }


   int DifferentialEvolution::mutatePopulation_Target_2_Bounded(double F1, double F2, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	    
	    int r3=rand() % (sizePopulation-3);
	    if(r3>=i) r3++;
	    if(r3>=r1) r3++;
	    if(r3>=r2) r3++;
	     

	    int r4=rand() % (sizePopulation-4);
	    if(r4>=i) r4++;
	    if(r4>=r1) r4++;
	    if(r4>=r2) r4++;
	    if(r4>=r3) r4++;

		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j])+F2*(Population[r3]->vectorCharacters[j] - Population[r4]->vectorCharacters[j]);

		if(nextPopulation[i]->vectorCharacters[j]<0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])<decoder->boundAttributes[2*j]){
			minF=((decoder->boundAttributes[2*j]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
			//cout<<minF<<"\t";
		}
		else if(nextPopulation[i]->vectorCharacters[j]>0 && (Population[i]->vectorCharacters[j] + minF * nextPopulation[i]->vectorCharacters[j])>decoder->boundAttributes[2*j+1]){
			minF=((decoder->boundAttributes[2*j+1]-Population[i]->vectorCharacters[j]) / nextPopulation[i]->vectorCharacters[j]);
			//minF=floor(minF*TRUNCCASE)/TRUNCCASE;
			minF=nexttoward(minF,0L);
		}
	
//if(minF<0)cout<<"a";	

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }


	}
	
	return 1;
      }
            
      int DifferentialEvolution::recombinePopulation(double CR){
	for(int i=0;i<sizePopulation;i++){
	    int jRand=rand() % nextPopulation[i]->sizeVec;
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){
	        int choice=rand() % 100;
	        
	        if(choice > CR && j != jRand)
			nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j];

	    }
	}
               
	return 1;
      }

      int DifferentialEvolution::recombinePopulationExp(){
	for(int i=0;i<sizePopulation;i++){
	    int init=rand() % (nextPopulation[i]->sizeVec-1);
	    int end=init + rand() % (nextPopulation[i]->sizeVec - init);
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){
	        
	        if( j < init || j > end)
			nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j];

	    }
	}
               
	return 1;
      }
      
int c=0;
      int DifferentialEvolution::selectPopulation(){
	Solution *swap;
	//cout<<"Iteracao  "<< c++<<"\n";
	for(int i=0;i<sizePopulation;i++){
	      decoder->decodifySolution(nextPopulation[i]);
	//cout<<Population[i]->vectorCharacters[0]<<"\t"<<Population[i]->vectorCharacters[1]<<"\t"<<Population[i]->score<<"\n";
	}
		//cout<<"\n----------------------------------------------------------\n";
	penalty->updatePenalty(Population,nextPopulation,sizePopulation,sizePopulation,decoder);
	for(int i=0;i<sizePopulation;i++){
	      if(!penalty->compareSolutions(Population[i],nextPopulation[i],decoder)){
		swap=Population[i];
		Population[i]=nextPopulation[i];
		nextPopulation[i]=swap;
			
			//cout<<*Population[i];
		if((best==NULL || penalty->compareSolutions(Population[i],best,decoder))){
		      if(!best) delete best;
		      best=Population[i]->clone();
		      UPLevelCallsBest=decoder->function->getUPLevelCalls();
		}
	      }
	}
               

	return 1;
      }
      
      int DifferentialEvolution::clearPopulation(){
	for(int i=0;i<sizePopulation; delete Population[i], delete nextPopulation[i++]);
	delete Population;
	delete nextPopulation;
	
	return 1;
      }
      

      
    
