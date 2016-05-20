#ifndef DIFFERENTIALEVOLUTION_INCLUDED
#define DIFFERENTIALEVOLUTION_INCLUDED 
#include "InputFunction.h"
#include "Solution.h"
#include "LagrangeMultp.h"
#include "LagrangeMultpSimplex.h"
#include "LagrangeMultpSmooth.h"
#include "APMPenalty.h"
#include "APMDEBPenalty.h"
#include "DEBPenalty.h"

class DifferentialEvolution{
      private: static Solution **Population;
      private: static Solution **nextPopulation;
      public: static Solution *best;
      public: static int UPLevelCallsBest;
      public: static int LWLevelSimplexCallsBest;
      private: static int sizePopulation;
      public: static SolutionDecoder *decoder;
      public: static PenaltySolution *penalty;

      public: static int initPopulation(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop);
      public: static int initPopulationNelderMeadMethod(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop); 

      
      public: static int mutatePopulation_Rand_1(double F, int begin, int end);
      public: static int mutatePopulation_Rand_2(double F, int begin, int end);
      public: static int mutatePopulation_RandToBest_1(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Rand_1_Bounded(double F, int begin, int end);      
      public: static int mutatePopulation_Rand_2_Bounded(double F, int begin, int end);      
      public: static int mutatePopulation_RandToBest_1_Bounded(double F1, double F2, int begin, int end);

      
      public: static int mutatePopulation_TargetToRand_1(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_TargetToBest_1(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Target_1(double F1, int begin, int end);
      public: static int mutatePopulation_Target_2(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_TargetToRand_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_TargetToBest_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Target_1_Bounded(double F1, int begin, int end);
      public: static int mutatePopulation_Target_2_Bounded(double F1, double F2, int begin, int end);

      
      public: static int mutatePopulation_BestToRand_1(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Best_1(double F1, int begin, int end);
      public: static int mutatePopulation_Best_2(double F1, double F2, int begin, int end);      
      public: static int mutatePopulation_BestToRand_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Best_1_Bounded(double F1, int begin, int end);
      public: static int mutatePopulation_Best_2_Bounded(double F1, double F2, int begin, int end);
      
      public: static int recombinePopulation(double CR, int begin, int end);
      public: static int recombinePopulationExp(double CR,int begin, int end);
      public: static int recombinePopulationSwap(int begin, int end);
      public: static int selectPopulation();     
      

      public: static int improveInitSet(int sizePopSearch, int sizePopNextStep, double find1, double find2);
      public: static int improveInitSetSimilarity(int sizePopSearch, int sizePopNextStep, double find1, double find2, int mutOption, double crossRate,int crossOpt, int sizeBest);
      
      public: static int improveInitSetDispersion(int sizePopSearch, int sizePopNextStep, double find1, double find2, int mutOption, double crossRate,int crossOpt, int sizeBest);
      
      
      public: static int clearPopulation();
      
      
      public: static int decodifyPopulation();
};


#endif
