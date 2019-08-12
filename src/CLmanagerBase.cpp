#include "CLmanagerBase.h"
#include <vector>
#include <exception>

using namespace clde;

std::list<CLmanagerBase::DeviceData> CLmanagerBase::getAllDevices(const cl_device_type& type)
{
    std::list<DeviceData> result;
    cl_uint number_of_platforms;
    cl_uint number_of_devices;
    cl_int err = 0;
    std::vector<cl_device_id> devices;
    clGetPlatformIDs(0, nullptr, &number_of_platforms);

    std::vector<cl_platform_id> platform_ids(number_of_platforms);
    clGetPlatformIDs(number_of_platforms, platform_ids.data(), nullptr);

    for(auto platform_it = platform_ids.begin(); platform_it != platform_ids.end(); ++platform_it)
    {
        err = clGetDeviceIDs(*platform_it, type, 0, nullptr, &number_of_devices);
        if (err < 0 )
            continue;
        devices.resize(number_of_devices);
        clGetDeviceIDs(*platform_it, type, number_of_devices, devices.data(), nullptr);
        for(auto device_it = devices.begin(); device_it != devices.end(); ++device_it)
        {
            result.push_back(DeviceData(*device_it, *platform_it));
        }

    }
    return result;
}

CLmanagerBase::~CLmanagerBase()
{
    if(m_queue)
        clReleaseCommandQueue(*m_queue);
    if(m_context)
        clReleaseContext(*m_context);
}

const cl_device_id& CLmanagerBase::device()
{
    if(!m_device_id)
    {
        auto devices = getAllDevices(CL_DEVICE_TYPE_GPU);
        if(devices.size() == 0)
        {
            devices = getAllDevices(CL_DEVICE_TYPE_CPU);
            if(devices.size() == 0)
                throw std::runtime_error("OpenCL: can't find any device");
        }
        m_device_id = devices.front();
    }
    return m_device_id->device;
}

const cl_context& CLmanagerBase::context()
{
    if(!m_context)
    {
        cl_int err = 0;
        if(!m_device_id)
            device();
        cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(m_device_id->platform), 0};
        m_context = clCreateContext(properties, 1, & (device()), nullptr, nullptr, &err);
        if(err < 0 )
            throw std::runtime_error("OpenCL: unable to create context");
    }
    return *m_context;
}

const cl_command_queue& CLmanagerBase::command_queue()
{
    if(!m_queue)
    {
        cl_int err = 0;
        if(!m_context)
            context();
        m_queue = clCreateCommandQueueWithProperties(*m_context, m_device_id->device, nullptr, &err);
        if(err < 0 )
            throw std::runtime_error("OpenCL: unable to create command queue");
    }
    return *m_queue;
}
