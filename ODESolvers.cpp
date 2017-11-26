// Implements the Runge-Kutta  ODE solver
// that derives from AbstractOdeSolver

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
RungeKuttaSolver::RungeKuttaSolver ( std::vector<double> (* f) (std::vector<double> initState, double t), double t0, double t1, std::vector<double> y0, double h, std::string filename) {

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
    
    // Store the coefficients used in the calculation (e.g. k)    
    double coeffs [n_eqs];

    // Sum of all k_n used to calculate the last step
    // (e.g k_0 +2k_1 + 2k_2 + k_3)
    double coeffSum [n_eqs];

    // Temporary state that is the parameter for the 
    // system  
    std::vector <double> tempState (n_eqs);

    // First step
    std::vector <double> ddt = SystemEquations(state, t);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];     
        tempState[i] = state[i] + 0.5*coeffs[i];
        coeffSum[i] = coeffs[i];
    }

    // Second
    ddt = SystemEquations(tempState, t+0.5*m_stepSize);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        tempState[i] += 0.5*coeffs[i];
        coeffSum[i] += coeffs[i]*2;
    }
   
    // Third
    ddt = SystemEquations(tempState, t+0.5*m_stepSize); 
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        tempState[i] += coeffs[i];
        coeffSum[i] += coeffs[i]*2;
    }

    // Fourth
    ddt = SystemEquations(tempState, t+m_stepSize);
    for (int i = 0; i < n_eqs; i++) {
        coeffs[i] = m_stepSize * ddt[i];
        coeffSum[i] += coeffs[i];
    }

    // Final
    for (int i = 0; i < n_eqs; i++) {
        state[i] += coeffSum[i]/6.0;

    }

}

