// Abstract class definition for an ODE solver

#include <vector>
#include <string>
#ifndef ABSTRACTODESOLVER_H
#define ABSTRACTODESOLVER_H

class AbstractOdeSolver
{
    protected:
        double stepSize;
        double initTime;
        double finalTime;
        
        std::string filename;       
        std::vector <double> initValues;
        std::vector <double> (*f) (std::vector <double>, double);
        virtual std::vector <double> SystemEquations(std::vector<double> state, double t);
        
        virtual void iterate (std::vector<double> &, double t)=0;
    public:
        void setStepSize(double h); 
        void setTimeInterval(double t0, double t1); 
        void setInitialValue(std::vector<double> initValues); 
        virtual double SolveEquation();
};

#endif
