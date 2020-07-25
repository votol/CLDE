#include "DERunge4.h"
#include <exception>
#include <iostream>

using namespace clde;


DERunge4::DERunge4(const std::shared_ptr<ICLmanager>& contex, IDEOperator* op):
    m_dt(1.0), m_Nsteps(1), m_Noutputs(1),m_blas(contex),
    m_operator(op), m_init(op->dimension().in_dim, contex),
    m_vec1(op->dimension().in_dim, contex),m_vec2(op->dimension().in_dim, contex),
    m_vec3(op->dimension().in_dim, contex),m_vec4(op->dimension().in_dim, contex),
    m_parameters(1)
{
    if(op->dimension().in_dim != op->dimension().out_dim)
        throw std::runtime_error("DERunge4 : input operator should be square");
}

void DERunge4::SetTimeStep(const double & in)
{
    if(in <= 0.0)
        m_dt = 1.0;
    else
        m_dt = in;
}

void DERunge4::SetStepsNumber(const unsigned int & in)
{
    m_Nsteps = in;
}

void DERunge4::SetOutputSteps(const unsigned int & in)
{
    m_Noutputs = in;
}

void DERunge4::SetParameters(const std::vector<double> & in)
{
    m_parameters = in;
}

void DERunge4::SetInitState(const std::vector<double> & in)
{
    if(in.size() != m_operator->dimension().in_dim)
        throw std::invalid_argument("Dimension of input vector should agree with operator");
    auto tmp = in;
    tmp[0] = 1.0;
    m_init = tmp;
}

void DERunge4::SetOutputs(const std::list<IDEOutput *> & in)
{
    m_outputs = in;
}

void DERunge4::calculate()
{
    double coe = double(m_Nsteps - 1) * m_dt/ double(m_Noutputs - 1);
    double output_threshold = 0.0;
    double koe;
    koe=0.0;
    m_parameters[0] = 0.0;
    m_vec1 = m_init;
    for(unsigned int ind = 0 ; ind < m_Nsteps; ind++)
    {
        if(m_parameters[0] >= output_threshold)
        {
            for(auto it = m_outputs.begin();it!=m_outputs.end();++it)
                (*it)->apply(m_vec1,m_parameters);
            output_threshold += coe;
        }

        m_vec2 = m_vec1;
        m_operator->apply(m_vec1, m_vec3, m_parameters);
        koe=m_dt/6.0;
        m_blas.Daxpy(m_vec2, m_vec3, koe);
        m_vec4 = m_vec1;
        koe=m_dt/2.0;
        m_blas.Daxpy(m_vec4, m_vec3, koe);
        m_parameters[0] += m_dt/2.0;
        m_operator->apply(m_vec4, m_vec3, m_parameters);
        koe=m_dt/3.0;
        m_blas.Daxpy(m_vec2, m_vec3, koe);
        m_vec4 = m_vec1;
        koe=m_dt/2.0;
        m_blas.Daxpy(m_vec4, m_vec3, koe);
        m_operator->apply(m_vec4, m_vec3, m_parameters);
        koe=m_dt/3.0;
        m_blas.Daxpy(m_vec2, m_vec3, koe);
        m_vec4 = m_vec1;
        koe=m_dt;
        m_blas.Daxpy(m_vec4, m_vec3, koe);
        m_parameters[0] += m_dt/2.0;
        m_operator->apply(m_vec4, m_vec3, m_parameters);
        koe=m_dt/6.0;
        m_blas.Daxpy(m_vec2, m_vec3, koe);
        m_vec1 = m_vec2;
    }
}
