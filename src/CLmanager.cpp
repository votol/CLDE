#include "CLmanager.h"
#include <string>

CLmanager::CLmanager(YAML::Node&& config)
{
    cl_device_type device_type;
    int device_num;
    if (YAML::Node type = config["cl_device_type"])
    {
        std::string tmp_device_type =  type.as<std::string>();

        if(tmp_device_type == "GPU")
            device_type = CL_DEVICE_TYPE_GPU;
        else if(tmp_device_type == "CPU")
            device_type = CL_DEVICE_TYPE_CPU;
        else
            device_type = CL_DEVICE_TYPE_GPU;
    }
    else
    {
        device_type = CL_DEVICE_TYPE_GPU;
    }

    if(YAML::Node number = config["cl_device_number"])
    {
        device_num = number.as<int>();
    }
    else
    {
        device_num = 0;
    }


    auto devices = getAllDevices(device_type);
    if(static_cast<size_t>(device_num) >= devices.size())
        throw std::runtime_error("OpenCL: can't find device");

    auto it = devices.begin();
    for (int counter = 0; counter < device_num; ++counter)
    {
        ++it;
    }
    m_device_id = *it;

    cl_int err;
    cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(m_device_id->platform), 0};
    m_context = clCreateContext(properties, 1, & (device()), nullptr, nullptr, &err);
    if(err < 0 )
        throw std::runtime_error("OpenCL: unable to create context");
    m_queue = clCreateCommandQueueWithProperties(*m_context, m_device_id->device, nullptr, &err);
    if(err < 0 )
        throw std::runtime_error("OpenCL: unable to create command queue");
}

const cl_device_id& CLmanager::device()
{
    return m_device_id->device;
}

const cl_context& CLmanager::context()
{
    return *m_context;
}

const cl_command_queue& CLmanager::command_queue()
{
    return *m_queue;
}

