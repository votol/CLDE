#pragma once
#include "ICLmanager.h"
#include "CLDataStorage.h"
#include <memory>

namespace clde
{
class CLBLAS
{
private:
    std::shared_ptr<ICLmanager> m_context;
    cl_kernel m_kernel_Daxpy;
public:
    CLBLAS(const std::shared_ptr<ICLmanager>& context);
    ~CLBLAS();
    //doing y = alpha*x + y for double variables
    void Daxpy(const CLDataStorage<double>& y, const CLDataStorage<double>& x, const double& alpha);
};
}
