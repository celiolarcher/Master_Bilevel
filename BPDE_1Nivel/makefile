flag=-O2

all: BPDE
BPDE: InputFunction.o Solution.o DifferentialEvolution.o LagrangeMultpAPM.o LagrangeMultpAPMSmooth.o LagrangeMultpFeasibilityComparison.o LagrangeMultpAPMSimplex.o main.o 
	g++ -o BPDE main.o Solution.o InputFunction.o LagrangeMultpAPM.o LagrangeMultpAPMSmooth.o LagrangeMultpFeasibilityComparison.o DifferentialEvolution.o LagrangeMultpAPMSimplex.o $(flag)

Solution.o: Solution.cpp Solution.h 
	g++ -c Solution.cpp $(flag)
	
InputFunction.o: InputFunction.cpp InputFunction.h ListFunction.h 
	g++ -c InputFunction.cpp $(flag)

DifferentialEvolution.o: DifferentialEvolution.cpp DifferentialEvolution.h
	g++ -c DifferentialEvolution.cpp $(flag)

main.o: main.cpp
	g++ -c main.cpp $(flag)

LagrangeMultpAPM.o: LagrangeMultpAPM.cpp LagrangeMultpAPM.h SolutionDecoder.h
	g++ -c LagrangeMultpAPM.cpp $(flag)
	
LagrangeMultpAPMSmooth.o: LagrangeMultpAPMSmooth.cpp LagrangeMultpAPMSmooth.h SolutionDecoder.h
	g++ -c LagrangeMultpAPMSmooth.cpp $(flag)
	
LagrangeMultpFeasibilityComparison.o: LagrangeMultpFeasibilityComparison.cpp LagrangeMultpFeasibilityComparison.h SolutionDecoder.h
	g++ -c LagrangeMultpFeasibilityComparison.cpp $(flag)

LagrangeMultpAPMSimplex.o: LagrangeMultpAPMSimplex.cpp LagrangeMultpAPMSimplex.h SolutionDecoder.h
	g++ -c LagrangeMultpAPMSimplex.cpp $(flag)
	
clean: 	
	rm *.o
	rm BPDE
