// Runge-Kutta method class 
// Other methods may be derived from AbstractOdeSolver

#include <vector>
#include "AbstractODESolver.h"

#ifndef ODESOLVERS_H
#define ODESOLVERS_H


typedef std::vector<double> sys_state; 
class RungeKuttaSolver : public AbstractOdeSolver {

    public:
        RungeKuttaSolver (std::vector<double> (* f) (std::vector<double>, double t),
                double t0, double t1, std::vector<double> initState, double h,
                 std::string filename = "");

    private:
        virtual void iterate (std::vector<double> &, double t);


};

#endif
