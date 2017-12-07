// Tests the HHSim class

#include<iostream>
#include<vector>
#include "ODESolvers.h"
#include <cmath>
#include "HHSim.h"

double currentStep(double t) {

    if (t >= 10 && t <= 60)
        return 15;

    return 0;

}

int main() {

    // Set initial conditions
    sys_state initial (4);
    initial[0] = 15;
    initial[1] = an(0)/(an(0) + bn(0));
    initial[2] = am(0)/(am(0) + bm(0));
    initial[3] = ah(0)/(ah(0) + bh(0));
    
    // Make new instance of HHSim 
    HHSim test(0.0, 20.0, initial, 0.005, "hh0005.txt");   
    test.SolveEquation(); 

    initial[0] = 0;
    HHSim testStep(0.0, 100.0, initial, 0.01, "hh-step.txt", currentStep);
    testStep.SolveEquation();
    return 0;

}
