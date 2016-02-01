#include "DifferentialEvolution.h" 
#include "stdlib.h"


 int validate(Solution *sol, InputFunction *function){
	if(function->getNEQConstraintNumberLW()>0){
		double list[function->getNEQConstraintNumberLW()];
		function->constraintsValueNEQLW(sol->vet,sol->vet+function->getDimensionUP(),list);
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		      if(list[i]>10e-5){
			  return 0;
		      }
		}
        }
        
        if(function->getEQConstraintNumberLW()>0){
		double list[function->getEQConstraintNumberLW()];
		function->constraintsValueEQLW(sol->vet,sol->vet+function->getDimensionUP(),list);
		for(int i=0;i<function->getEQConstraintNumberLW();i++){
		      if(list[i]>10e-5){
			  return 0;
		      }
		}
        }
        
        
        if(function->getNEQConstraintNumberUP()>0){
		double list[function->getNEQConstraintNumberUP()];
		function->constraintsValueNEQUP(sol->vet,sol->vet+function->getDimensionUP(),list);
		for(int i=0;i<function->getNEQConstraintNumberUP();i++){
		      if(list[i]>10e-5){
			  return 0;
		      }
		}
        }
        
        
        if(function->getEQConstraintNumberUP()>0){
		double list[function->getEQConstraintNumberUP()];
		function->constraintsValueEQUP(sol->vet,sol->vet+function->getDimensionUP(),list);
		for(int i=0;i<function->getEQConstraintNumberUP();i++){
		      if(list[i]>10e-5){
			  return 0;
		      }
		}
        }

        if(function->getKKTConstraintNumber()>0){
		double list[function->getKKTConstraintNumber()];
		function->constraintsValueKKT(sol->vet,sol->vet+function->getDimensionUP(),sol->vet+function->getDimensionUP()+function->getDimensionLW(),sol->vet+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
		for(int i=0;i<function->getKKTConstraintNumber();i++){
		      if(list[i]<-10e-8 || list[i]>10e-8){
			  return 0;
		      }
		}
        }

        if(function->getNEQConstraintNumberLW()>0){
		double list[function->getNEQConstraintNumberLW()];
		function->constraintsSlackness(sol->vet,sol->vet+function->getDimensionUP(),sol->vet+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
		for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		      if(list[i]<-10e-8 || list[i]>10e-8){
			  return 0;
		      }
		}
        }

	return 1;
     }

























      Solution **DifferentialEvolution::Population;
      Solution **DifferentialEvolution::nextPopulation;
      Solution *DifferentialEvolution::best=NULL;
      int DifferentialEvolution::sizePopulation;


      int DifferentialEvolution::initPopulation(InputFunction *function, int sizePop){
	Population=new Solution*[sizePop];
	nextPopulation=new Solution*[sizePop];
	sizePopulation=sizePop;
	
	for(int i=0;i<sizePopulation; Population[i]=new Solution(function), nextPopulation[i]=new Solution(function), Population[i++]->initRandom(function));  //alterar initRandom com sizeVet e bounds
        
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
	     
	    
	    for(int j=0;j<nextPopulation[i]->sizeVet;j++){    
	        nextPopulation[i]->vet[j]=Population[r1]->vet[j] + F*(Population[r2]->vet[j] - Population[r3]->vet[j]);
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
	    int jRand=rand() % nextPopulation[i]->sizeVet;
	    
	    for(int j=0;j<nextPopulation[i]->sizeVet;j++){
	        int choice=rand() % 100;
	        
	        if(choice > CR && j != jRand)
		nextPopulation[i]->vet[j]=Population[i]->vet[j];
	    }
	}
               
	return 1;
      }
      
#include <iostream>
using namespace std;
      int DifferentialEvolution::selectPopulation(InputFunction *function){
	Solution *swap;
	
	for(int i=0;i<sizePopulation;i++){
	      int val=validate(Population[i],function);
	      int valNext=validate(nextPopulation[i],function);
	      if((!val && valNext) || ((val && valNext || !val) && Population[i]->calcScore(function) > nextPopulation[i]->calcScore(function))){//avaliacao(Population[i]) > avaliacao(nextPopulation[i])){    //determinar avaliação
		swap=Population[i];
		Population[i]=nextPopulation[i];
		nextPopulation[i]=swap;
	      }
	      if(valNext || val && ( best==NULL || Population[i]->calcScore(function) < best->calcScore(function))) best=Population[i];
	}
               
	return 1;
      }
      
      int DifferentialEvolution::clearPopulation(){
	for(int i=0;i<sizePopulation; delete Population[i], delete nextPopulation[i++]);
	delete Population;
	delete nextPopulation;
	
	return 1;
      }





    
