flag=-O2

all: BPDE
BPDE: InputFunction.o Solution.o DifferentialEvolution.o OptLangrageMultp.o main.o 
	g++ -o BPDE main.o Solution.o InputFunction.o OptLangrageMultp.o DifferentialEvolution.o $(flag)

Solution.o: Solution.cpp Solution.h 
	g++ -c Solution.cpp $(flag)
	
InputFunction.o: InputFunction.cpp InputFunction.h ListFunction.h 
	g++ -c InputFunction.cpp $(flag)

DifferentialEvolution.o: DifferentialEvolution.cpp DifferentialEvolution.h
	g++ -c DifferentialEvolution.cpp $(flag)

main.o: main.cpp
	g++ -c main.cpp $(flag)

OptLangrageMultp.o: OptLangrageMultp.cpp OptLangrageMultp.h
	g++ -c OptLangrageMultp.cpp $(flag)

clean: 	
	rm *.o
	rm BPDE
