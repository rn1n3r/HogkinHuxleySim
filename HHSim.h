// Class extending RungeKuttaSolver to represent
// the Hodgkin Huxley model specifically

#include "ODESolvers.h"
#include <vector>

#ifndef HHSIM_H
#define HHSIM_H

double an(double v);
double am(double v);
double ah(double v);
double bn(double v);
double bm(double v);
double bh(double v);

class HHSim : public RungeKuttaSolver {


    public:
        HHSim ( double t0, double t1, std::vector<double> initState, double h,
                 std::string filename = "", double (*I)(double) = 0, double gNa = 120, double gK = 36,
             double gL = 0.3, double VNa = 115, double VK = -12, 
             double VL = 10.613, double C = 1);

    private:
        std::vector<double> SystemEquations (std::vector <double>, double t);
        
        // Conductances
        double gNa; double gK; double gL;
        double VNa; double VK; double VL;

        double C;
        double (*I) (double); 

};

#endif
