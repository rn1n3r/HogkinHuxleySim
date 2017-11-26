// Abstract class for an ODE solver
// Functions to set step size, time interval
// and initial values

#include "AbstractODESolver.h"
#include <fstream>
#include <iostream>
#include <vector>
// Setter function for step size
void AbstractOdeSolver::setStepSize(double h) { 
    stepSize = h; 
}

// Setter function for initial and final time point
void AbstractOdeSolver::setTimeInterval(double t0, double t1) {
    initTime = t0; 
    finalTime = t1;
};

// Setter function for the initial value
void AbstractOdeSolver::setInitialValue(std::vector<double> initValues ) {
    AbstractOdeSolver::initValues = initValues;
}

std::vector<double> AbstractOdeSolver::SystemEquations (std::vector<double> state, double t) {

    return AbstractOdeSolver::f(state, t);
}

// Solve the ODE
// Calculates the steps for the interval and iterates through
// the time interval. Outputs the solution to a textfile 
// specified by the member function "filename"

// Since the iteration function is implemented as a member virtual
// function, the SolveEquation method can be used by any 
// inherited class that might use a different method
double AbstractOdeSolver::SolveEquation () {
    double t = initTime;
    std::vector<double> state = initValues;

    double sqerror = 0;

    std::ofstream output(filename);

    // Calculate the number of steps 
    int steps =  (int)((finalTime - initTime)/stepSize);

    // Use a for loop over the integer value for number
    // of steps to avoid accumulating error (in t)
    for (int i = 0; i <= steps; i++) {

        if (!filename.empty()) {
            output << t << " ";
            for (int j = 0; j < state.size(); j++) {
                output << state[j];
                if (j < state.size() - 1)
                    output << " ";
            }
            output << std::endl;
        }
        else {
            std::cout << t << " ";
            for (int j = 0; j < state.size(); j++) {
                std::cout  << state[j];
                if (j < state.size() - 1)
                    std::cout << " ";
            }
            std::cout << std::endl;
        }
        // Iterate 
        iterate(state, t);
        t += stepSize;
    }

    // Output message
    if (!filename.empty()) {
        std::cout << "Solution written to " << filename << std::endl;
    }
    output.close();
    return 0;
}
