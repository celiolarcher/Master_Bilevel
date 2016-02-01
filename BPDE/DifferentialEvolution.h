#ifndef DIFFERENTIALEVOLUTION_INCLUDED
#define DIFFERENTIALEVOLUTION_INCLUDED 
#include "InputFunction.h"
#include "Solution.h"

class DifferentialEvolution{
      private: static Solution **Population;
      private: static Solution **nextPopulation;
      public: static Solution *best;
      private: static int sizePopulation;

      public: static int initPopulation(InputFunction *function, int sizePop); 
      public: static int mutatePopulation(double F);
      public: static int recombinePopulation(double CR);
      public: static int selectPopulation(InputFunction *function);     
      
      public: static int clearPopulation();
};


#endif
