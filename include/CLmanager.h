#pragma once
#include "CLmanagerBase.h"
#include "yaml-cpp/yaml.h"

class CLmanager: public CLmanagerBase
{
public:
    CLmanager(YAML::Node&&);
    virtual ~CLmanager() override = default;
    virtual const cl_device_id& device() override;
    virtual const cl_context& context() override;
    virtual const cl_command_queue& command_queue() override;
};
