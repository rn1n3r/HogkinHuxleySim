#include<cmath>
#include "HHSim.h"
#include <vector>

double am (double v) {
    return 0.1 * (25 - v) / (exp(0.1 * (25 - v)) - 1); 
}

double bm (double v) {
    return 4 * exp(-v / 18);
}

double ah (double v) {
    return 0.07 * exp(-v / 20);
}

double bh (double v) {
    return 1 / (exp((30 - v) / 10) + 1);
}

double an (double v) {
    return 0.01 * (-v + 10) / (exp(0.1*(-v + 10)) - 1);
}

double bn (double v) {
    return 0.125 * exp(-v/80);
}

double noCurrent (double t) {
    return 0;
}

std::vector<double> HHSim::SystemEquations(std::vector<double> vnmh, double t) {

    const double I = 0; 

    sys_state ddt(4);
    ddt[0] = (-gK * pow(vnmh[1], 4) * (vnmh[0] - VK) 
             -gNa * pow(vnmh[2], 3) * vnmh[3]* (vnmh[0] - VNa)
             -gL *(vnmh[0] - VL) + I)/C;

    ddt[1] = an(vnmh[0]) * (1 - vnmh[1]) - bn(vnmh[0]) * vnmh[1];
    ddt[2] = am(vnmh[0]) * (1 - vnmh[2]) - bm(vnmh[0]) * vnmh[2];
    ddt[3] = ah(vnmh[0]) * (1 - vnmh[3]) - bh(vnmh[0]) * vnmh[3];
    return ddt;
}


HHSim::HHSim(double t0, double t1, std::vector<double> initState, double h,
             std::string filename, double gNa, double gK,
             double gL, double VNa, double VK, 
             double VL, double C) : RungeKuttaSolver(NULL, t0, t1, initState, h, filename) 

{    
    HHSim::gNa = gNa;
    HHSim::gK = gK;
    HHSim::gL = gL;
    HHSim::VNa = VNa;
    HHSim::VL = VL;
    HHSim::VK = VK;
    HHSim::C = C;
}

