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
      int DifferentialEvolution::LWLevelSimplexCallsBest;
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
		LWLevelSimplexCallsBest=decoder->function->getLWLevelSimplexCalls();
	    }
	}
	
	return 1;
      }
      

      
      void quickSortSolution(Solution *eliteSet[], int left, int right, SolutionDecoder *decoder){  //Utilizado para ordenar as soluções por melhor score.
          int i, j;
          Solution *y;
          Solution *x;
          i = left;
          j = right;
          x = eliteSet[(left + right) / 2];
          while(i <= j){
	  while(DifferentialEvolution::penalty->compareSolutions(eliteSet[i],x,decoder) && i < right)  i++;
	  while(DifferentialEvolution::penalty->compareSolutions(x,eliteSet[j],decoder) && j > left)   j--;
	  if(i <= j){
	      y = eliteSet[i];
	      eliteSet[i] = eliteSet[j];
	      eliteSet[j] = y;
	      i++;
	      j--;
	  }
          }
          if(j > left)  quickSortSolution(eliteSet, left, j,decoder);
          if(i < right) quickSortSolution(eliteSet,  i, right,decoder);
      }

      
      
      int DifferentialEvolution::initPopulationNelderMeadMethod(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop){
	double alpha=1, gama=2, theta=0.5;
        
        
	Population=new Solution*[sizePop];
	nextPopulation=new Solution*[sizePop];
	sizePopulation=sizePop;
	DifferentialEvolution::decoder=decoder;
	DifferentialEvolution::penalty=penalty;
	
	Solution *initPop[sizePop];

	
	for(int i=0;i<sizePopulation; initPop[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), nextPopulation[i]=new Solution(decoder->solutionSize,decoder->constraintsNumber), initPop[i++]->initRandom(decoder->boundAttributes)); 
        
	
	for(int i=0;i<sizePopulation;i++){
	    decoder->decodifySolution(initPop[i]);
	}
	
	penalty->updatePenalty(initPop,NULL,sizePopulation,0,decoder);

	
	Solution *simplexGroup[sizePop];
	
	for(int i=0;i<sizePop;i++){
	      
	      int mark[sizePop];
	      for(int j=0;j<sizePop;mark[j++]=0);
	      
	      Solution *pointsSelectd[decoder->solutionSize+1];
	      int worst=0;
	      int best=0;
	      
	      for(int j=0;j<decoder->solutionSize+1;j++){
	        
		 int r=rand() % (sizePopulation-j);
		 
		 int count=0,k=0;
		 for(k=0;k<sizePopulation &&  count<r;k++){
		    if(!mark[k])count++;
		}
		
		pointsSelectd[j]=initPop[k];
		
		if(penalty->compareSolutions(pointsSelectd[worst],pointsSelectd[j],decoder)) worst=j;
		
		if(penalty->compareSolutions(pointsSelectd[j],pointsSelectd[best],decoder)) best=j;
	      }
	    
	      
	      double centroidPoint[decoder->solutionSize];
	      
	      for(int j=0;j<decoder->solutionSize;j++){
		centroidPoint[j]=0;
		for(int k=0;k<decoder->solutionSize+1;k++){
		      if(k!=worst)
		          centroidPoint[j]+=pointsSelectd[j]->vectorCharacters[k];
		}
		centroidPoint[j]/=decoder->solutionSize;
	      }
	      
	      double xAux[decoder->solutionSize];
	      
    	      for(int j=0;j<decoder->solutionSize;j++){
	          xAux[j]=centroidPoint[j] + alpha*(centroidPoint[j]-pointsSelectd[worst]->vectorCharacters[j]);
	      }

	      Solution *reflect=new Solution(decoder->solutionSize,decoder->constraintsNumber);
	      reflect->initValue(xAux);
	      
	      decoder->decodifySolution(reflect);
	      
	      if(penalty->compareSolutions(reflect,pointsSelectd[best],decoder)){
		  for(int j=0;j<decoder->solutionSize;j++){
		      xAux[j]=centroidPoint[j] + gama*(reflect->vectorCharacters[j] - centroidPoint[j]);
		  }
		  Solution *expand=new Solution(decoder->solutionSize,decoder->constraintsNumber);
		  expand->initValue(xAux);
		  
		  decoder->decodifySolution(expand);
		  
		  if(penalty->compareSolutions(expand,pointsSelectd[best],decoder)){
		  //if(penalty->compareSolutions(expand,reflect,decoder)){

		      simplexGroup[i]=expand;
		      delete reflect;
		  }
		  else{
		      simplexGroup[i]=reflect;
		      delete expand;
		  }
	      }else  if(penalty->compareSolutions(reflect,pointsSelectd[worst],decoder)){
	        
		  for(int j=0;j<decoder->solutionSize;j++){
		      xAux[j]=centroidPoint[j] + theta*(pointsSelectd[worst]->vectorCharacters[j] - centroidPoint[j]);
		  }
		  
		  Solution *contract=new Solution(decoder->solutionSize,decoder->constraintsNumber);
		  contract->initValue(xAux);
		  
		  decoder->decodifySolution(contract);
		  
		  if(penalty->compareSolutions(contract,pointsSelectd[worst],decoder)){
		  //if(penalty->compareSolutions(contract,reflect,decoder)){
		      simplexGroup[i]=contract;
		      delete reflect;
		  }
		  else{
		      simplexGroup[i]=reflect;
		      delete contract;
		  }
	      }else{  //Criterio?
	        	
		 double rd= ((double) rand()/ RAND_MAX) * (1+gama+alpha);
		 
		 for(int j=0;j<decoder->solutionSize;j++){
		        
		        //xAux[j]=centroidPoint[j] + theta*(pointsSelectd[worst]->vectorCharacters[j] - centroidPoint[j]);
		        xAux[j]=pointsSelectd[worst]->vectorCharacters[j] + rd*(centroidPoint[j] - pointsSelectd[worst]->vectorCharacters[j]);
		 }
		 
		 Solution *randPoint=new Solution(decoder->solutionSize,decoder->constraintsNumber);
		 randPoint->initValue(xAux);
		 
		 decoder->decodifySolution(randPoint);
		 
		 simplexGroup[i]=randPoint;
	      }
	      
	}
	
	
	penalty->updatePenalty(initPop,simplexGroup,sizePopulation,sizePopulation,decoder);
		
	quickSortSolution(initPop,0,sizePopulation-1,decoder);
	quickSortSolution(simplexGroup,0,sizePopulation-1,decoder);
	
	int countInitPop=0, countSimplex=0;
	for(int i=0; i<sizePopulation ; i++){
	    if(penalty->compareSolutions(initPop[countInitPop],simplexGroup[countSimplex],decoder)){
	        Population[i]=initPop[countInitPop++];
	    }else{
	        Population[i]=simplexGroup[countSimplex++];
	    }
	}
	
	for(int i=countInitPop;i<sizePopulation;i++) delete initPop[i];
	for(int i=countSimplex;i<sizePopulation;i++) delete simplexGroup[i];
	

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

		for(int i=0;i<10e5 && decoder->function->getLWLevelSimplexCalls()<5*10e6 && feasibleSolutions<sizePopNextStep/4; i++){
			DifferentialEvolution::mutatePopulation_Target_2_Bounded(find1,find2,0,sizePopSearch-feasibleSolutions);
			DifferentialEvolution::recombinePopulationSwap(0,sizePopulation);
			DifferentialEvolution::selectPopulation();

			for(int j=0;j<sizePopSearch-feasibleSolutions && feasibleSolutions<sizePopNextStep/4;j++){
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
	
		for(int j=sizePopNextStep-feasibleSolutions;j<sizePopulation;j++) delete Population[j];

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
      extern bool InfeasibleAvaliation;

      
      int DifferentialEvolution::improveInitSetSimilarity(int sizePopSearch, int sizePopNextStep, double find1, double find2, int mutOption, double crossRate,int crossOpt, int sizeBest){
		int MIN_POP=6;
        
		Solution **populationNew=new Solution *[sizePopNextStep+sizeBest];
		int feasibleSolutions=0;		

		for(int i=0;i<10e5 && decoder->function->getLWLevelSimplexCalls()<5*10e6/2 && feasibleSolutions<sizePopNextStep && sizePopulation>=MIN_POP; i++){
			double rd= ((double) rand()/ RAND_MAX) * (find2);
		  
			/*
			if(mutOption==1)
			    DifferentialEvolution::mutatePopulation_Rand_1_Bounded(find1,0,sizePopulation);
			else if(mutOption==2)
			    DifferentialEvolution::mutatePopulation_Target_2_Bounded(find1,find2,0,sizePopulation);
			else if(mutOption==3)
			    DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(find1,find2,0,sizePopulation);
			else
			    DifferentialEvolution::mutatePopulation_TargetToBest_1_Bounded(find1,find2,0,sizePopulation);
			*/
						
			if(mutOption==1)
			    DifferentialEvolution::mutatePopulation_Rand_1_Bounded(find1+rd,0,sizePopulation);
			else if(mutOption==2)
			    DifferentialEvolution::mutatePopulation_Target_2_Bounded(find1+rd,find1+rd,0,sizePopulation);
			else if(mutOption==3)
			    DifferentialEvolution::mutatePopulation_TargetToRand_1_Bounded(find1+rd,find1+rd,0,sizePopulation);
			
			
			
			if(crossOpt==1)
			    DifferentialEvolution::recombinePopulationSwap(0,sizePopulation);
			else if(crossOpt==2)
			    DifferentialEvolution::recombinePopulation(crossRate,0,sizePopulation);	   
			else if(crossOpt==4)
			    DifferentialEvolution::recombinePopulationExp(crossRate,0,sizePopulation);	   
			
			
			DifferentialEvolution::selectPopulation();

			for(int j=0;j<sizePopulation && feasibleSolutions<sizePopNextStep;j++){
				if(Population[j]->feasible){
					populationNew[feasibleSolutions]=Population[j];
					Population[j]=Population[sizePopulation-1];
					sizePopulation--;
					feasibleSolutions++;
					
					for(int del=0;del<sizePopSearch/sizePopNextStep-1 && sizePopulation>MIN_POP;del++){
					    double minDiff=10e5;
					    int deleted=0;
					    for(int k=0;k<sizePopulation;k++){

					          //double aux=populationNew[feasibleSolutions-1]->diffSquareSolution(Population[k]);

					          double aux=populationNew[feasibleSolutions-1]->diffMaxSolution(Population[k]);
					      
					          if(aux<minDiff){
						  minDiff=aux;
						  deleted=k;
					          }
					    }
					    delete Population[deleted];
					    Population[deleted]=Population[sizePopulation-1];
					    sizePopulation--;
					}
					
				}
			}
		}

		for(int j=feasibleSolutions;j<sizePopNextStep;j++) populationNew[j]=Population[j-feasibleSolutions];
	
		for(int j=sizePopNextStep-feasibleSolutions;j<sizePopulation;j++) delete Population[j];

		Solution **nextPopulationNew=new Solution*[sizePopNextStep+sizeBest];
		
		
		
		
		
		//for(int j=sizePopNextStep;j<sizePopNextStep+sizeBest;j++) populationNew[j]=best->clone();  //Melhorar??
		
		for(int j=sizePopNextStep;j<sizePopNextStep+sizeBest; populationNew[j]=new Solution(decoder->solutionSize,decoder->constraintsNumber), populationNew[j]->initRandom(decoder->boundAttributes),decoder->decodifySolution(populationNew[j++])); 
		

		
		
		
		
		for(int j=0;j<sizePopNextStep+sizeBest && j<sizePopSearch;j++)nextPopulationNew[j]=nextPopulation[j];
		
		for(int j=sizePopSearch;j<sizePopNextStep+sizeBest;j++)  nextPopulationNew[j]=new Solution(decoder->solutionSize,decoder->constraintsNumber);


		for(int j=sizePopNextStep+sizeBest;j<sizePopSearch;j++) delete nextPopulation[j];

		delete Population;
		delete nextPopulation;
		
		Population=populationNew;
		nextPopulation=nextPopulationNew;
		sizePopulation=sizePopNextStep+sizeBest;
	/*	
		
		delete penalty;
		
		penalty = new APMPenalty();

		InfeasibleAvaliation=1;

		
		for(int i=0;i<sizePopulation;i++) decoder->decodifySolution(Population[i]);
	*/
	
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
      
      
      int DifferentialEvolution::mutatePopulation_Rand_2(double F,int begin, int end){
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
	    
	    int r5=rand() % (sizePopulation-5);
	    if(r5>=i) r5++;
	    if(r5>=r1) r5++;
	    if(r5>=r2) r5++;
	    if(r5>=r3) r5++;
	    if(r5>=r4) r5++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=F;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	            nextPopulation[i]->vectorCharacters[j]=(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]) + (Population[r4]->vectorCharacters[j] - Population[r5]->vectorCharacters[j]);
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	//	if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) cout<<Population[r1]->vectorCharacters[j]<<"\t"<<minF<<"\t"<<(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j])<<"\t"<<nextPopulation[i]->vectorCharacters[j]<<"\n";
	//	if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) cout<<Population[r1]->vectorCharacters[j]<<"\t"<<minF<<"\t"<<(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j])<<"\t"<<nextPopulation[i]->vectorCharacters[j]<<"\n";

//		if(nextPopulation[i]->vectorCharacters[j]<decoder->boundAttributes[2*j]) cout<<nextPopulation[i]->vectorCharacters[j]-decoder->boundAttributes[2*j]<<"\n";
//		if(nextPopulation[i]->vectorCharacters[j]>decoder->boundAttributes[2*j+1]) cout<<nextPopulation[i]->vectorCharacters[j]-decoder->boundAttributes[2*j]<<"\n";
	        

	    }
	    
	    //      if(nextPopulation[i]->diffSquareSolution(Population[i])<10e-5) cout<<"a";

//if(minF<0)cout<<*nextPopulation[i];	

	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_TargetToRand_1(double F1, double F2, int begin, int end){
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
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }


	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_TargetToBest_1(double F1, double F2, int begin, int end){
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

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }


      int DifferentialEvolution::mutatePopulation_BestToRand_1(double F1, double F2, int begin, int end){
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

	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=best->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }

      
       int DifferentialEvolution::mutatePopulation_Best_1(double F1, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]);
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=best->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_Best_2(double F1, double F2, int begin, int end){
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
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]) + F2*(Population[r3]->vectorCharacters[j] - Population[r4]->vectorCharacters[j]);
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=best->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];

	    }
	}
	
	return 1;
      }
      
      
      
      int DifferentialEvolution::mutatePopulation_Target_1(double F1, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;

		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	          nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]);
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }


	}
	
	return 1;
      }
      

      int DifferentialEvolution::mutatePopulation_Target_2(double F1, double F2, int begin, int end){
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
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]) + F2*(Population[r3]->vectorCharacters[j] - Population[r4]->vectorCharacters[j]);
	    }
	    
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }

	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_RandToBest_1(double F1, double F2,int begin, int end){
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
	    }
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=Population[r1]->vectorCharacters[j] + minF*nextPopulation[i]->vectorCharacters[j];
	    }

	}
	
	return 1;
      }
      
      
      int DifferentialEvolution::mutatePopulation_Rand_1_Bounded(double F,int begin, int end){
	for(int i=begin;i<end;i++){    //int flag;
	//  do{flag=0;
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
	  //  	        if(Population[r1]->diffMaxSolution(nextPopulation[i])<10e-5)flag=1;

//	          if(nextPopulation[i]->diffSquareSolution(Population[i])<10e-5) cout<<"a";
	//  }while(flag);

	}
	
	return 1;
      }
      
      
      
      int DifferentialEvolution::mutatePopulation_Rand_2_Bounded(double F,int begin, int end){
	for(int i=begin;i<end;i++){ //   int flag;
	 // do{flag=0;
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
	    
	    int r5=rand() % (sizePopulation-5);
	    if(r5>=i) r5++;
	    if(r5>=r1) r5++;
	    if(r5>=r2) r5++;
	    if(r5>=r3) r5++;
	    if(r5>=r4) r5++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=F;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	            nextPopulation[i]->vectorCharacters[j]=(Population[r2]->vectorCharacters[j] - Population[r3]->vectorCharacters[j]) + (Population[r4]->vectorCharacters[j] - Population[r5]->vectorCharacters[j]);

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
	   // 	        if(Population[r1]->diffMaxSolution(nextPopulation[i])<10e-3)flag=1;

	    
//	          if(nextPopulation[i]->diffSquareSolution(Population[i])<10e-5) cout<<"a";
	 // }while(flag);
	    
	    //      if(nextPopulation[i]->diffSquareSolution(Population[i])<10e-5) cout<<"a";

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

      
       int DifferentialEvolution::mutatePopulation_Best_1_Bounded(double F1, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;
	     
		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]);

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
      
      
      int DifferentialEvolution::mutatePopulation_Best_2_Bounded(double F1, double F2, int begin, int end){
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
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]) + F2*(Population[r3]->vectorCharacters[j] - Population[r4]->vectorCharacters[j]);

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
      
      
      
      int DifferentialEvolution::mutatePopulation_Target_1_Bounded(double F1, int begin, int end){
	for(int i=begin;i<end;i++){    
	    int r1=rand() % (sizePopulation-1);
	    if(r1>=i) r1++;
	    
	    int r2=rand() % (sizePopulation-2);
	    if(r2>=i) r2++;
	    if(r2>=r1) r2++;

		//cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";

	    double minF=1;

	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){    
	          nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]);

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
	        nextPopulation[i]->vectorCharacters[j]=F1*(Population[r1]->vectorCharacters[j] - Population[r2]->vectorCharacters[j]) + F2*(Population[r3]->vectorCharacters[j] - Population[r4]->vectorCharacters[j]);

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
            
      int DifferentialEvolution::recombinePopulation(double CR ,int begin, int end){
	CR*=100;
	for(int i=begin;i<end;i++){
	    int jRand=rand() % nextPopulation[i]->sizeVec;
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){
	        int choice=rand() % 100;
	        
	        if(choice > CR && j != jRand)
		    nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j];

	    }
	}
               
	return 1;
      }

      int DifferentialEvolution::recombinePopulationSwap(int begin, int end){
	for(int i=begin;i<end;i++){
	    int init=rand() % (nextPopulation[i]->sizeVec);
	    int end=init + rand() % (nextPopulation[i]->sizeVec - init);
	    
	    for(int j=0;j<nextPopulation[i]->sizeVec;j++){
	        
	        if( j < init || j > end)
			nextPopulation[i]->vectorCharacters[j]=Population[i]->vectorCharacters[j];

	    }
	}
               
	return 1;
      }
      
      
      int DifferentialEvolution::recombinePopulationExp(double CR, int begin, int end){
	CR*=100;
	for(int i=begin;i<end;i++){
	    int init=rand() % (nextPopulation[i]->sizeVec);
	    int j=init;
	    int L=0,choice;

	    
	    /*
	    double aux[nextPopulation[i]->sizeVec];
	    
	    for(int k=0;k<nextPopulation[i]->sizeVec;k++) aux[k]=Population[i]->vectorCharacters[k];
	    
	    do{
	        aux[j]=nextPopulation[i]->vectorCharacters[j];
	        j=(j+1) % (nextPopulation[i]->sizeVec);
	        L++;
	        choice=rand() % 100;
	    }while(choice<CR && L<nextPopulation[i]->sizeVec);
	    
	    
	    for(int k=0;k<nextPopulation[i]->sizeVec;k++) nextPopulation[i]->vectorCharacters[k]=aux[k];
*/
	    
	    
	    do{
	        j=(j+1) % (nextPopulation[i]->sizeVec);
	        L++;
	        choice=rand() % 100;
	    }while(choice<CR && L<nextPopulation[i]->sizeVec);
	    
	    if(j>init){
	        for(int k=0;k<init;k++)
	          nextPopulation[i]->vectorCharacters[k]=Population[i]->vectorCharacters[k];
	        for(int k=j;k<nextPopulation[i]->sizeVec;k++)		
	          nextPopulation[i]->vectorCharacters[k]=Population[i]->vectorCharacters[k];
	    }else{
	        for(int k=j;k<init;k++)
	          nextPopulation[i]->vectorCharacters[k]=Population[i]->vectorCharacters[k];
	    }
	    
	    
	}
               
	return 1;
      }
      
      
int c=0;
      int DifferentialEvolution::selectPopulation(){
	Solution *swap;
	//cout<<"Iteracao  "<< c++<<"\n";
/*
	for(int i=0;i<sizePopulation;i++){
	  int flag=0;
	  for(int j=0;j<sizePopulation;j++){
	      if(nextPopulation[i]->diffMaxSolution(Population[j])<10e-5){
	          flag=1;
	          delete nextPopulation[i];
	          nextPopulation[i]=Population[j]->clone();
	      }
	  }
	
	  if(!flag)	decoder->decodifySolution(nextPopulation[i]);
	  
	//cout<<Population[i]->vectorCharacters[0]<<"\t"<<Population[i]->vectorCharacters[1]<<"\t"<<Population[i]->score<<"\n";
	}
	
	*/
	
	for(int i=0;i<sizePopulation;i++){
	    if(nextPopulation[i]->diffMaxSolution(Population[i])<10e-5){
	      delete nextPopulation[i];
	      nextPopulation[i]=Population[i]->clone();
	    }else   decoder->decodifySolution(nextPopulation[i]);
	    //for(int j=0;j<Population[i]->sizeVec;j++)
	      //  cout<<Population[i]->vectorCharacters[j]<<"\t";
	     //cout<<"\t"<<Population[i]->score<<"\n";
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
		//      cout<<decoder->function->getUPLevelCalls()<<"\t"<<best->score<<"\t"<<Population[i]->score<<"\t"<<best->diffSquareSolution(Population[i])<<"\n";
		      if(!best) delete best;
		      best=Population[i]->clone();
		      UPLevelCallsBest=decoder->function->getUPLevelCalls();
		      LWLevelSimplexCallsBest=decoder->function->getLWLevelSimplexCalls();
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
      
         
      int DifferentialEvolution::decodifyPopulation(){
	  for(int i=0;i<sizePopulation;i++) decoder->decodifySolution(Population[i]);
	  decoder->decodifySolution(best);
	  return 1;
      }

      
    
