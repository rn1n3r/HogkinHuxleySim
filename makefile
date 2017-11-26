all: bin/main

bin/main: AbstractODESolver.cpp ODESolvers.cpp main.cpp 
	g++ -std=c++11 AbstractODESolver.cpp ODESolvers.cpp main.cpp -o bin/main




