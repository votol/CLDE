#pragma once

#include "IFuncCalculator.h"

namespace clde
{
class FakeFuncCalculator: public IFuncCalculator
{
private:
    CLDataStorage<double> m_data;
    unsigned int m_ind = 0;
public:
    FakeFuncCalculator(const std::shared_ptr<ICLmanager>& context);
    virtual void init(const size_t& num) override;
    virtual const CLDataStorage<double>& data() const override;
    virtual const unsigned int& process(const double&) override;
};
}
