#pragma once
#include "CLDataStorage.h"
#include "CLBLAS.h"
#include "IDEOperator.h"
#include "IDEOutput.h"
#include <list>

namespace clde{
class DERunge4
{
private:
    double m_dt;
    unsigned int m_Nsteps;
    unsigned int m_Noutputs;

    CLBLAS m_blas;
    IDEOperator* m_operator;
    std::list<IDEOutput*> m_outputs;

    CLDataStorage<double> m_init;
    CLDataStorage<double> m_vec1;
    CLDataStorage<double> m_vec2;
    CLDataStorage<double> m_vec3;
    CLDataStorage<double> m_vec4;


    std::vector<double> m_parameters;
public:
    DERunge4(const std::shared_ptr<ICLmanager>&, IDEOperator*);
    void SetTimeStep(const double&);
    void SetStepsNumber(const unsigned int&);
    void SetOutputSteps(const unsigned int&);
    void SetParameters(const std::vector<double>&);
    void SetInitState(const std::vector<double>&);
    void SetOutputs(const std::list<IDEOutput*> &);
    void calculate();
};

};
