#include "DifferentialEvolution.h" 
#include "stdlib.h"


#include <iostream>
using namespace std;
      Solution **DifferentialEvolution::Population;
      Solution **DifferentialEvolution::nextPopulation;
      Solution *DifferentialEvolution::best=NULL;
      int DifferentialEvolution::sizePopulation;
      SolutionDecoder *DifferentialEvolution::decoder;


      int DifferentialEvolution::initPopulation(SolutionDecoder *decoder, int sizePop){
	Population=new Solution*[sizePop];
	nextPopulation=new Solution*[sizePop];
	sizePopulation=sizePop;
	DifferentialEvolution::decoder=decoder;
	
	for(int i=0;i<sizePopulation; Population[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), nextPopulation[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), Population[i++]->initRandom(decoder->boundAttributes)); 
        
	for(int i=0;i<sizePopulation;i++){
	    decoder->decodifySolution(Population[i]);
		//cout<<*Population[i];
	    if(Population[i]->feasible && (best==NULL || decoder->compareSolutions(Population[i],best))){
	        if(!best) delete best;
	        best=Population[i]->clone();
	    }
	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation(double F){
	for(int i=0;i<sizePopulation;i++){    
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
/*
		if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) nextPopulation[i]->vectorCharacters[j]=decoder->boundAttributes[2*j];
		if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) nextPopulation[i]->vectorCharacters[j]=decoder->boundAttributes[2*j+1];
*/
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
      
int c=0;
      int DifferentialEvolution::selectPopulation(){
	Solution *swap;
	//cout<<"Iteracao  "<< c++<<"\n";
	for(int i=0;i<sizePopulation;i++){
	      decoder->decodifySolution(nextPopulation[i]);
	//cout<<Population[i]->vectorCharacters[0]<<"\t"<<Population[i]->vectorCharacters[1]<<"\t"<<Population[i]->score<<"\n";
	}
	//cout<<"\n----------------------------------------------------------\n";
	decoder->updatePenalty(Population,nextPopulation,sizePopulation,sizePopulation);
	for(int i=0;i<sizePopulation;i++){
	      if(!decoder->compareSolutions(Population[i],nextPopulation[i])){
		swap=Population[i];
		Population[i]=nextPopulation[i];
		nextPopulation[i]=swap;
			
		//cout<<*Population[i];
		if(Population[i]->feasible && (best==NULL || decoder->compareSolutions(Population[i],best))){
		      if(!best) delete best;
		      best=Population[i]->clone();
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
      

      
    
