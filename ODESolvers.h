// Runge-Kutta method class 
// Other methods may be derived from AbstractOdeSolver

#include <vector>
#include "AbstractODESolver.h"


class RungeKuttaSolver : public AbstractOdeSolver {

    public:
        RungeKuttaSolver (std::vector<double> (* SystemEquations) (std::vector<double>, double t),
                double t0, double t1, std::vector<double> initState, double h,
                 std::string filename = "");

    private:
        virtual void iterate (std::vector<double> &, double t);


};
