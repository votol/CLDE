#pragma once
#include <optional>
#include <list>
#include "ICLmanager.h"


class CLmanagerBase: public ICLmanager
{
protected:
    struct DeviceData
    {
        cl_device_id device;
        cl_platform_id platform;
        DeviceData(const cl_device_id& d, const cl_platform_id& p):device(d), platform(p){}
    };
    std::optional<DeviceData> m_device_id = std::nullopt;
    std::optional<cl_context> m_context = std::nullopt;
    std::optional<cl_command_queue> m_queue = std::nullopt;
    std::list<DeviceData> getAllDevices(const cl_device_type&);

public:
    CLmanagerBase() = default;
    virtual ~CLmanagerBase() override;
    virtual const cl_device_id& device() override;
    virtual const cl_context& context() override;
    virtual const cl_command_queue& command_queue() override;
};
