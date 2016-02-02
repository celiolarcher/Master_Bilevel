#include "DifferentialEvolution.h" 
#include "stdlib.h"
#include "OptLangrageMultp.h"


      Solution **DifferentialEvolution::Population;
      Solution **DifferentialEvolution::nextPopulation;
      Solution *DifferentialEvolution::best=NULL;
      int DifferentialEvolution::sizePopulation;


      int DifferentialEvolution::initPopulation(InputFunction *function, int sizePop){
	Population=new Solution*[sizePop];
	nextPopulation=new Solution*[sizePop];
	sizePopulation=sizePop;
	int sizeVector=function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW();
	int sizeConstraints=function->getEQConstraintNumberUP()+function->getNEQConstraintNumberUP()+function->getEQConstraintNumberLW()+function->getNEQConstraintNumberLW()+function->getKKTConstraintNumber()+function->getNEQConstraintNumberLW();
	
	double bounds[2*sizeVector];
	for(int i=0;i<function->getDimensionUP()+function->getDimensionLW();i++){
	    bounds[2*i]=function->bounds[2*i];
	    bounds[2*i+1]=function->bounds[2*i+1];
	}
	
	for(int i=function->getDimensionUP()+function->getDimensionLW();i<sizeVector;i++){
	    bounds[2*i]=0;
	    bounds[2*i+1]=20;
	}
	
	for(int i=0;i<sizePopulation; Population[i]=new Solution(sizeVector,sizeConstraints), nextPopulation[i]=new Solution(sizeVector,sizeConstraints), Population[i++]->initRandom(bounds)); 
        
	for(int i=0;i<sizePopulation;i++){
	    validateSolution(Population[i],function);
	    if(Population[i]->feasible && (best==NULL || compareSolutions(Population[i],best,function)))
	        best=Population[i];
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
	     
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + F*(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]);
	    }
	
	/*	int rad=rand()%100;
	    if(rad<10){
			int rand2=rand()%4;
			nextPopulation[i]->vet[2+rand2]=0;
		}
*/
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
      
#include <iostream>
using namespace std;
      int DifferentialEvolution::selectPopulation(InputFunction *function){
	Solution *swap;
	
	updateAPM(Population,sizePopulation,function);
	
	for(int i=0;i<sizePopulation;i++){
	      validateSolution(nextPopulation[i],function);
	      if(!compareAPM(Population[i],nextPopulation[i],function)){
		swap=Population[i];
		Population[i]=nextPopulation[i];
		nextPopulation[i]=swap;
			
		if(Population[i]->feasible && (best==NULL || compareSolutions(Population[i],best,function)))
		      best=Population[i];
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
      

      
    
