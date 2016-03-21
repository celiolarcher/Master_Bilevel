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
      private: static int sizePopulation;
      public: static SolutionDecoder *decoder;
      public: static PenaltySolution *penalty;

      public: static int initPopulation(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop); 
      public: static int mutatePopulation_Rand_1(double F, int begin, int end);
      public: static int mutatePopulation_Rand_1_Bounded(double F, int begin, int end);
      public: static int mutatePopulation_RandToBest_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_TargetToRand_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_TargetToBest_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_BestToRand_1_Bounded(double F1, double F2, int begin, int end);
      public: static int mutatePopulation_Target_2_Bounded(double F1, double F2, int begin, int end);
      public: static int recombinePopulation(double CR);
      public: static int recombinePopulationExp();
      public: static int selectPopulation();     
      

      public: static int improveInitSet(int sizePopSearch, int sizePopNextStep, double find1, double find2);
      public: static int clearPopulation();
};


#endif
