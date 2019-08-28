#pragma once

#include "IFuncCalculator.h"

namespace clde
{
class WienerFuncCalculator: public IFuncCalculator
{
private:
    std::shared_ptr<ICLmanager> m_context;
    CLDataStorage<double> m_data;
    CLDataStorage<double> m_random_numbers;
    std::optional<cl_kernel> m_kernel;
    size_t m_threads_num =0;
    size_t m_work_group = 32;
    unsigned int m_random_shift = 0;
    unsigned int m_ind = 0;
    double m_current_time = 0.0;
    double m_delta_time = 0.0;
public:
    WienerFuncCalculator(const std::shared_ptr<ICLmanager>& context);
    virtual ~WienerFuncCalculator() override;
    virtual void init(const size_t& num) override;
    virtual const CLDataStorage<double>& data() const override;
    virtual const unsigned int& process(const double&) override;
};
}
