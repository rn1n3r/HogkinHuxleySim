#include<iostream>
#include<vector>
#include "ODESolvers.h"
#include <cmath>

typedef std::vector<double> sys_state;

sys_state lorenz (sys_state xyz, double t) {
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0/3.0;
    sys_state ddt(3); 
    ddt[0] = sigma*(xyz[1] - xyz[0]);
    ddt[1] = R * xyz[0] - xyz[1] - xyz[0] * xyz[2];
    ddt[2] = xyz[0]*xyz[1] - b*xyz[2];
    return ddt;
}

sys_state f(sys_state xyz, double t) {

    sys_state ddt(2);
    ddt[0] = xyz[1];
    ddt[1] = xyz[0]*6.0 - xyz[1];
        
    return ddt;
}
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

sys_state hh(sys_state vnmh, double t) {

    const double gNa = 120;
    const double gK = 36;
    const double gL = 0.3;
    const double VNa = 115;
    const double VK = -12;
    const double VL = 10.613;
    const double C = 1;
    double I = 0;
    

    sys_state ddt(4);
    ddt[0] = (-gK * pow(vnmh[1], 4) * (vnmh[0] - VK) 
             -gNa * pow(vnmh[2], 3) * vnmh[3]* (vnmh[0] - VNa)
             -gL *(vnmh[0] - VL) + I)/C;

    ddt[1] = an(vnmh[0]) * (1 - vnmh[1]) - bn(vnmh[0]) * vnmh[1];
    ddt[2] = am(vnmh[0]) * (1 - vnmh[2]) - bm(vnmh[0]) * vnmh[2];
    ddt[3] = ah(vnmh[0]) * (1 - vnmh[3]) - bh(vnmh[0]) * vnmh[3];
    return ddt;
}

int main() {

    sys_state initial (4);
    initial[0] = 20.0;
    initial[1] = an(0)/(an(0) + bn(0));
    initial[2] = am(0)/(am(0) + bm(0));
    initial[3] = ah(0)/(ah(0) + bh(0));
    
    RungeKuttaSolver test (hh, 0.0, 20.0, initial, 0.01, "hh.txt");   
    
    sys_state initialLorenz (3);
    initialLorenz[0] = initialLorenz[1] = initialLorenz[2] = 10.0;
    RungeKuttaSolver lorenzTest (lorenz, 0.0, 10.0, initialLorenz, 0.01, "lorenz.txt"); 
    //sys_state initial (2);
    //initial[0] = 3;
    //initial[1] = 1; 
    //RungeKuttaSolver test (f, 0, 1, initial, 0.01);
    test.SolveEquation(); 
    lorenzTest.SolveEquation();
    return 0;


}
