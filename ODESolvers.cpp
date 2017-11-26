// Assignment 4, Problem 2
// Implements the forward Euler ODE solver
// and the Runge-Kutta method by defining
// classes that inherit from AbstractOdeSolver

#include "AbstractODESolver.h"
#include "ODESolvers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <string>
#include <vector>

// Constructor
// f = specified function representing the ODE
// t0 and t1 = time interval
// y0 = initial value y(t=0)
// filename = name of text file to output results to 
// t = specified true solution. 
RungeKuttaSolver::RungeKuttaSolver ( std::vector<double> (* f) (std::vector<double> initState, double t), double t0,
        double t1, std::vector<double> y0, double h, std::string filename) {

    // Store the specified function  
    RungeKuttaSolver::SystemEquations = f;

    // Set time interval, step, and initial value
    setStepSize (h);
    setTimeInterval (t0, t1);
    setInitialValue (y0);

    // For writing to file
    RungeKuttaSolver::filename = filename;

}



// Overload of the iterate function
void RungeKuttaSolver::iterate (std::vector<double> &state, double t) {
    int n_eqs = state.size();
    
    // First step
    double coeffs [n_eqs];
    double coeffSum [n_eqs]; 
    std::vector <double> state2 (n_eqs);
    std::vector <double> ddt = SystemEquations(state, t);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];     
        state2[i] = state[i] + 0.5*coeffs[i];
        coeffSum[i] = coeffs[i];
    }
    

    // Second
    ddt = SystemEquations(state2, t+0.5*m_stepSize);
    std::vector <double> state3 (n_eqs);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        state3[i] = state2[i] + 0.5*coeffs[i];
        coeffSum[i] += coeffs[i]*2;
    }
   
    // Third
    ddt = SystemEquations(state3, t+0.5*m_stepSize); 
    std::vector <double> state4 (n_eqs);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        state4[i] = state3[i] + coeffs[i];
        coeffSum[i] += coeffs[i]*2;
    }

    // Fourth
    ddt = SystemEquations(state4, t+m_stepSize);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        coeffSum[i] += coeffs[i];
    }

    // Final
    for (int i = 0; i < n_eqs; i++) {
        state[i] += coeffSum[i]/6.0;

    }

}





