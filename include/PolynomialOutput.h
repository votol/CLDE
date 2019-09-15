#pragma once
#include "IDEOutput.h"
#include "PolynomialOperator.h"

namespace clde
{
class PolynomialOutput: public IDEOutput
{
private:
    PolynomialOperator m_oper;
    std::shared_ptr<ICLmanager> m_context;
    CLDataStorage<double> m_vector;
    unsigned int m_num;
    unsigned int m_out_dim;
    unsigned int m_calc_count;

protected:
    std::vector<double> m_result;
    unsigned int m_current_output = 0;

public:
    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    PolynomialOutput(_InputIterator begin,  _InputIterator end,
                     const unsigned int& num, const OperatorDimension& dim, const unsigned int& calc_count,
                     const std::shared_ptr<ICLmanager>& context): m_oper(begin, end, num, dim, context), m_context(context),
                     m_vector(dim.out_dim * num + 1, context), m_num(num), m_out_dim(dim.out_dim), m_calc_count(calc_count), m_result(calc_count * m_out_dim, 0.0)
    {

    }
    ~PolynomialOutput() override;
    virtual void apply(const CLDataStorage<double>& , const std::vector<double>&) override;
};
};
