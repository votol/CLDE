#pragma once
#include <cstddef>
#include "CLDataStorage.h"

namespace clde
{

class IFuncCalculator{
public:
    virtual ~IFuncCalculator() = default;
    virtual void init(const size_t& num) = 0;
    virtual const CLDataStorage<double>& data() const = 0;
    virtual const unsigned int& process(const double& time) = 0;
};
}
