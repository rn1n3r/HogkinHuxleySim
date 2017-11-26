// Abstract class definition for an ODE solver

#include <vector>
#include <string>
#ifndef ABSTRACTODESOLVER_H
#define ABSTRACTODESOLVER_H

class AbstractOdeSolver
{
    protected:
        double m_stepSize;
        double m_initTime;
        double m_finalTime;
        
        std::string filename;       
        std::vector <double> m_initValues;
        std::vector <double> (*SystemEquations)(std::vector<double>, double);

        virtual void iterate (std::vector<double> &, double t)=0;
    public:
        void setStepSize(double h); 
        void setTimeInterval(double t0, double t1); 
        void setInitialValue(std::vector<double> initValues); 
        virtual double SolveEquation();
};

#endif
