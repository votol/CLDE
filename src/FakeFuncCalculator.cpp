#include "FakeFuncCalculator.h"
#include <vector>

using namespace clde;

FakeFuncCalculator::FakeFuncCalculator(const std::shared_ptr<ICLmanager>& context):
    m_data(context)
{

}

void FakeFuncCalculator::init(const size_t &num)
{
    m_data = std::vector<double>(num, 0.0);
}

const CLDataStorage<double>& FakeFuncCalculator::data() const
{
    return m_data;
}

const unsigned int& FakeFuncCalculator::process(const double &)
{
    return m_ind;
}
