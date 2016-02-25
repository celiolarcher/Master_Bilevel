#ifndef DIFFERENTIALEVOLUTION_INCLUDED
#define DIFFERENTIALEVOLUTION_INCLUDED 
#include "InputFunction.h"
#include "Solution.h"
#include "LagrangeMultpAPM.h"
#include "LagrangeMultpAPMSmooth.h"

class GeneticAlgorithmSimilarity{
      private: static Solution **Population;
      private: static Solution **nextPopulation;
      private: static Solution **eliteSet;
      public: static Solution *best;
      private: static int sizePopulation;
      public: static SolutionDecoder *decoder;

      public: static int initPopulation(SolutionDecoder *decoder, int sizePop); 
      public: static int mutatePopulation(double F);
      public: static int recombinePopulation(double CR);
      public: static int selectPopulation();     
      
      public: static int clearPopulation();
};


#endif
