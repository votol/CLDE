#include "PolynomialOutput.h"
#include <iostream>

using namespace clde;

void PolynomialOutput::apply(const CLDataStorage<double>& in, const std::vector<double>& param)
{
    if( m_current_output == m_calc_count)
        return;

    m_oper.apply(in, m_vector, param);
    auto result = m_vector.read();
    unsigned int count = 0;
    for(auto it = ++(result.begin()); it != result.end(); ++ it)
    {
        m_result[m_out_dim * m_current_output + count % m_out_dim] += *it;
        count++;
    }

    for(unsigned int ind = m_out_dim * m_current_output; ind < m_out_dim * m_current_output + m_out_dim; ++ ind)
    {
        m_result[ind] /= double(m_num);
    }
    m_current_output ++;
}

PolynomialOutput::~PolynomialOutput()
{

}
