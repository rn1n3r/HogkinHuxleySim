#include<iostream>
#include<vector>
#include "ODESolvers.h"
#include <cmath>
#include "HHSim.h"

//typedef std::vector<double> sys_state;

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

int main() {

    sys_state initial (4);
    initial[0] = 20.0;
    initial[1] = an(0)/(an(0) + bn(0));
    initial[2] = am(0)/(am(0) + bm(0));
    initial[3] = ah(0)/(ah(0) + bh(0));
    
    HHSim test(0.0, 20.0, initial, 0.01, "hh.txt");    
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
