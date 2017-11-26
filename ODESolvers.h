// Assignment 4, Problem 2
// Implements the forward Euler ODE solver
// and the Runge-Kutta method by defining
// classes that inherit from AbstractOdeSolver

#include <vector>
#include "AbstractODESolver.h"

// Runge-Kutta method class 
// Derived from ForwardEulerSolver since the 
// only difference is the iterate step

class RungeKuttaSolver : public AbstractOdeSolver {

    public:
        RungeKuttaSolver (std::vector<double> (* SystemEquations) (std::vector<double>, double t), double t0,
                double t1, std::vector<double> initState, double h, std::string filename = "");


    private:
        virtual void iterate (std::vector<double> &, double t);
        


        

};
