// Abstract class for an ODE solver
// Functions to set step size, time interval
// and initial values

#include "AbstractODESolver.h"
#include <fstream>
#include <iostream>

// Setter function for step size
void AbstractOdeSolver::setStepSize(double h) { 
    m_stepSize = h; 
}

// Setter function for initial and final time point
void AbstractOdeSolver::setTimeInterval(double t0, double t1) {
    m_initTime = t0; 
    m_finalTime = t1;
};

// Setter function for the initial value
void AbstractOdeSolver::setInitialValue(std::vector<double> initValues ) {
    m_initValues = initValues;
}

// Solve the ODE
// Calculates the steps for the interval and iterates through
// the time interval. Outputs the solution to a textfile 
// specified by the member function "filename"

// Since the iteration function is implemented as a member virtual
// function, the SolveEquation method can be used by any 
// inherited class that might use a different method
double AbstractOdeSolver::SolveEquation () {
    double t = m_initTime;
    std::vector<double> state = m_initValues;

    double sqerror = 0;

    std::ofstream output(filename);

    // Calculate the number of steps 
    int steps =  (int)((m_finalTime - m_initTime)/m_stepSize);

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
        t += m_stepSize;
    }

    // Output message
    if (!filename.empty()) {
        std::cout << "Solution written to " << filename << std::endl;
    }
    output.close();
    return 0;
}
