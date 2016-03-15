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
      private: static int sizePopulation;
      public: static SolutionDecoder *decoder;
      public: static PenaltySolution *penalty;

      public: static int initPopulation(SolutionDecoder *decoder, PenaltySolution *penalty,int sizePop); 
      public: static int mutatePopulation(double F);
      public: static int mutatePopulationBounded(double F);
      public: static int mutatePopulationBestBounded(double F1, double F2);
      public: static int mutatePopulationTargetBestBounded(double F1, double F2);
      public: static int mutatePopulationTargetBounded(double F1, double F2);
      public: static int recombinePopulation(double CR);
      public: static int selectPopulation();     
      
      public: static int clearPopulation();
};


#endif
