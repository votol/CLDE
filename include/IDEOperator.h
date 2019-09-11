#pragma once
#include <cstddef>
#include <vector>
#include "CLDataStorage.h"

namespace clde {
struct OperatorDimension
{
    size_t in_dim = 0;
    size_t out_dim = 0;
};

class IDEOperator
{
public:
    virtual ~IDEOperator() = default;
    virtual const OperatorDimension& dimension() = 0;
    virtual void apply(const CLDataStorage<double>& in, CLDataStorage<double>& out,
                                            const std::vector<double>&) = 0;
};
};
