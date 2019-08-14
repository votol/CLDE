#pragma once
#include <vector>
#include "CLDataStorage.h"

namespace clde{
class IDEOutput
{
public:
    virtual ~IDEOutput() = default;
    virtual void apply(const CLDataStorage<double>& in, const std::vector<double>&) = 0;
};

};
