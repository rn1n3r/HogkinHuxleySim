all: bin/main

bin/main: AbstractODESolver.cpp ODESolvers.cpp HHSim.cpp main.cpp 
	g++ -std=c++11 AbstractODESolver.cpp ODESolvers.cpp HHSim.cpp main.cpp -o bin/main




