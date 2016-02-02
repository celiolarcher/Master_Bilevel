#ifndef PENALTYFUNCTION_INCLUDED
#define PENALTYFUNCTION_INCLUDED  


    double validateFunction(double eqConstraintList[], int sizeEqConstraint, double neqConstraintList[], int sizeNeqConstraint,){
          if(function->getNEQConstraintNumberLW()>0){
		  double list[function->getNEQConstraintNumberLW()];
		  function->constraintsValueNEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),list);
		  for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		        if(list[i]>10e-5){
			    return 0;
		        }
		  }
          }
          
          if(function->getEQConstraintNumberLW()>0){
		  double list[function->getEQConstraintNumberLW()];
		  function->constraintsValueEQLW(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),list);
		  for(int i=0;i<function->getEQConstraintNumberLW();i++){
		        if(list[i]>10e-5){
			    return 0;
		        }
		  }
          }
          
          
          if(function->getNEQConstraintNumberUP()>0){
		  double list[function->getNEQConstraintNumberUP()];
		  function->constraintsValueNEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),list);
		  for(int i=0;i<function->getNEQConstraintNumberUP();i++){
		        if(list[i]>10e-5){
			    return 0;
		        }
		  }
          }
          
          
          if(function->getEQConstraintNumberUP()>0){
		  double list[function->getEQConstraintNumberUP()];
		  function->constraintsValueEQUP(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),list);
		  for(int i=0;i<function->getEQConstraintNumberUP();i++){
		        if(list[i]>10e-5){
			    return 0;
		        }
		  }
          }

          if(function->getKKTConstraintNumber()>0){
		  double list[function->getKKTConstraintNumber()];
		  function->constraintsValueKKT(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
		  for(int i=0;i<function->getKKTConstraintNumber();i++){
		        if(list[i]<-10e-8 || list[i]>10e-8){
			    return 0;
		        }
		  }
          }

          if(function->getNEQConstraintNumberLW()>0){
		  double list[function->getNEQConstraintNumberLW()];
		  function->constraintsSlackness(sol->vectorCharacters,sol->vectorCharacters+function->getDimensionUP(),sol->vectorCharacters+function->getDimensionUP()+function->getDimensionLW()+function->getEQConstraintNumberLW(),list);
		  for(int i=0;i<function->getNEQConstraintNumberLW();i++){
		        if(list[i]<-10e-8 || list[i]>10e-8){
			    return 0;
		        }
		  }
          }

          return 1;
     }











#endif
