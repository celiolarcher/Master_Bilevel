#ifndef DIFFERENTIALEVOLUTION_INCLUDED
#define DIFFERENTIALEVOLUTION_INCLUDED 
#include "InputFunction.h"
#include "Solution.h"
#include "LagrangeMultp.h"
#include "LagrangeMultpSimplex.h"
#include "LagrangeMultpSmooth.h"
#include "APMPenalty.h"
#include "APMDEBPenalty.h"
#include "LemkeLW.h"
#include "LemkeLWFillSolution.h"
//#include "DELowerLevel.h"
#include "DEBPenalty.h"


class DifferentialEvolution{
      public: Solution **Population;
      public: Solution **nextPopulation;
      public: Solution *best;
      public: int UPLevelCallsBest;
      public: int LWLevelSimplexCallsBest;
      private: int sizePopulation;
      public: SolutionDecoder *decoder;
      public: PenaltySolution *penalty;
      public: int bestSizePop;
      
      public: DifferentialEvolution(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop, int bestSizePop);

      public: int initPopulationNelderMeadMethod(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop); 

      
      public: int mutatePopulation_Rand_1(double F, int popElm);
      public: int mutatePopulation_Rand_2(double F, int popElm);
      public: int mutatePopulation_RandToBest_1(double F1, double F2, int popElm);
      public: int mutatePopulation_Rand_1_Wall(double F, int popElm);
      public: int mutatePopulation_Rand_2_Wall(double F, int popElm);
      public: int mutatePopulation_RandToBest_1_Wall(double F1, double F2, int popElm);
      public: int mutatePopulation_Rand_1_Bounded(double F, int popElm);      
      public: int mutatePopulation_Rand_2_Bounded(double F, int popElm);      
      public: int mutatePopulation_RandToBest_1_Bounded(double F1, double F2, int popElm);

      
      public: int mutatePopulation_TargetToRand_1(double F1, double F2, int popElm);
      public: int mutatePopulation_TargetToBest_1(double F1, double F2, int popElm);
      public: int mutatePopulation_Target_1(double F1, int popElm);
      public: int mutatePopulation_Target_2(double F1, double F2, int popElm);
      public: int mutatePopulation_TargetToRand_1_Wall(double F1, double F2, int popElm);
      public: int mutatePopulation_TargetToBest_1_Wall(double F1, double F2, int popElm);
      public: int mutatePopulation_Target_1_Wall(double F1, int popElm);
      public: int mutatePopulation_Target_2_Wall(double F1, double F2, int popElm);
      public: int mutatePopulation_TargetToRand_1_Bounded(double F1, double F2, int popElm);
      public: int mutatePopulation_TargetToBest_1_Bounded(double F1, double F2, int popElm);
      public: int mutatePopulation_Target_1_Bounded(double F1, int popElm);
      public: int mutatePopulation_Target_2_Bounded(double F1, double F2, int popElm);

      
      
      
      public: int mutatePopulation_BestToRand_1(double F1, double F2, int popElm);
      public: int mutatePopulation_Best_1(double F1, int popElm);
      public: int mutatePopulation_Best_2(double F1, double F2, int popElm);      
      public: int mutatePopulation_BestToRand_1_Wall(double F1, double F2, int popElm);
      public: int mutatePopulation_Best_1_Wall(double F1, int popElm);
      public: int mutatePopulation_Best_2_Wall(double F1, double F2, int popElm);      
      public: int mutatePopulation_BestToRand_1_Bounded(double F1, double F2, int popElm);
      public: int mutatePopulation_Best_1_Bounded(double F1, int popElm);
      public: int mutatePopulation_Best_2_Bounded(double F1, double F2, int popElm);
      
      public: int recombinePopulation(double CR, int popElm);
      public: int recombinePopulationExp(double CR,int popElm);
      public: int recombinePopulationSwap(int popElm);
      public: int selectPopulation(int popElm);     
      public: int selectPopulationBestPath(int popElm);

      public: int improveInitSetDE(int sizePopSearch, int sizePopNextStep, double find1, double find2);
      public: int improveInitSetSimilarityDE(int sizePopSearch, int sizePopNextStep, double find1, double find2, int mutOption, double crossRate,int crossOpt, int sizeBest, int intervalAddDE);
      
      public: int improveInitSetDispersionDE(int sizePopSearch, int sizePopNextStep, double find1, double find2, int mutOption, double crossRate,int crossOpt, int sizeBest, int intervalAddDE);

      public: int improveInitSetSimilarity(int sizePopSearch, int sizePopNextStep);
      
      
      
      public: int resetPopulation();
      public: int clearPopulation();
      public: int averagePopulation(double average[]);
      public: int findClosePopulation(double close[],Solution *sol);
      

      
      
      public: int decodifyPopulation();
      
      public: int calcVariance(double variance[]);

      
      
      public: int printPopulation();

};


#endif
