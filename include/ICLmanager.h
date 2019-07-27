#pragma once

#include <CL/cl.h>

class ICLmanager
{
public:
    virtual ~ICLmanager() = default;
    virtual const cl_device_id& device() = 0;
    virtual const cl_context& context() = 0;
};
