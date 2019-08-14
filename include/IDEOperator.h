#pragma once
#include <cstddef>
#include <vector>
#include "CLDataStorage.h"

namespace clde {
class IDEOperator
{
public:
    virtual ~IDEOperator() = default;
    virtual const size_t& dimension() = 0;
    virtual void apply(const CLDataStorage<double>& in, CLDataStorage<double>& out,
                                            const std::vector<double>&) = 0;
};
};
