flag=-O2

all: BPDE
BPDE: InputFunction.o Solution.o DifferentialEvolution.o  LagrangeMultp.o LagrangeMultpSimplex.o LagrangeMultpSmooth.o APMPenalty.o APMDEBPenalty.o DEBPenalty.o main.o 
	g++ -o BPDE main.o Solution.o InputFunction.o LagrangeMultp.o LagrangeMultpSimplex.o LagrangeMultpSmooth.o DifferentialEvolution.o APMPenalty.o APMDEBPenalty.o DEBPenalty.o  $(flag)

Solution.o: Solution.cpp Solution.h 
	g++ -c Solution.cpp $(flag)
	
InputFunction.o: InputFunction.cpp InputFunction.h ListFunction.h ListFunctionSMD.h ListFunctionJaq.h
	g++ -c InputFunction.cpp $(flag)

DifferentialEvolution.o: DifferentialEvolution.cpp DifferentialEvolution.h
	g++ -c DifferentialEvolution.cpp $(flag)

main.o: main.cpp DefaultParameters.h
	g++ -c main.cpp $(flag)

LagrangeMultp.o: LagrangeMultp.cpp LagrangeMultp.h SolutionDecoder.h
	g++ -c LagrangeMultp.cpp $(flag)

LagrangeMultpSimplex.o: LagrangeMultpSimplex.cpp LagrangeMultpSimplex.h SolutionDecoder.h
	g++ -c LagrangeMultpSimplex.cpp $(flag)

LagrangeMultpSmooth.o: LagrangeMultpSmooth.cpp LagrangeMultpSmooth.h SolutionDecoder.h
	g++ -c LagrangeMultpSmooth.cpp $(flag)

APMPenalty.o: APMPenalty.cpp APMPenalty.h Penalty.h
	g++ -c APMPenalty.cpp $(flag)

APMDEBPenalty.o: APMDEBPenalty.cpp APMDEBPenalty.h Penalty.h
	g++ -c APMDEBPenalty.cpp $(flag)

DEBPenalty.o: DEBPenalty.cpp DEBPenalty.h Penalty.h
	g++ -c DEBPenalty.cpp $(flag)

clean: 	
	rm *.o
	rm BPDE
