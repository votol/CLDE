#include <string>
#include <string_view>
#include <exception>
#include "CLBLAS.h"
#include "kernels/Daxpy.h"
using namespace clde;

cl_kernel createKernel(const char* source,  const std::string& func_name,const std::shared_ptr<ICLmanager>& context)
{
    cl_int err = 0;

    auto program_size = std::string_view(source).size();
    auto cl_program = clCreateProgramWithSource(context->context(), 1,
          &(source), &program_size, &err);

    if(err < 0)
    {
        throw std::runtime_error("OpenCL: unable to create program");
    }

    err = clBuildProgram(cl_program, 0, nullptr, nullptr, nullptr, nullptr);

    if(err < 0)
    {
        size_t log_size;
        clGetProgramBuildInfo(cl_program, context->device(), CL_PROGRAM_BUILD_LOG,
                    0, nullptr, &log_size);
        std::string program_log(log_size, ' ');
        clGetProgramBuildInfo(cl_program, context->device(), CL_PROGRAM_BUILD_LOG,
               log_size + 1, program_log.data(), nullptr);
        throw std::runtime_error("OpenCL: unable to build programm: " + program_log);
    }

    auto cl_kernel = clCreateKernel(cl_program, func_name.c_str(), &err);
    if(err < 0) {
        throw std::runtime_error("OpenCL: couldn't create a kernel");
    }

    clReleaseProgram(cl_program);
    return cl_kernel;
}


CLBLAS::CLBLAS(const std::shared_ptr<ICLmanager>& context):m_context(context)
{
    m_kernel_Daxpy = createKernel(Daxpy_kernel, "Daxpy", m_context);
}

CLBLAS::~CLBLAS()
{
    clReleaseKernel(m_kernel_Daxpy);
}

void CLBLAS::Daxpy(const CLDataStorage<double>& y, const CLDataStorage<double>& x, const double& alpha)
{
    if(y.size() != x.size())
        throw std::runtime_error("Daxpy: input vectors should have the same size");

    size_t num = y.size();

    cl_int err = 0;
    cl_event run_event;

    err = clSetKernelArg(m_kernel_Daxpy, 0, sizeof(cl_mem), &y.data());
    err |= clSetKernelArg(m_kernel_Daxpy, 1, sizeof(cl_mem), &x.data());
    err |= clSetKernelArg(m_kernel_Daxpy, 2, sizeof(double), &alpha);

    if(err < 0) {
        throw std::runtime_error("OpenCL: couldn't set an argument for the Daxpy kernel");
    }

    err = clEnqueueNDRangeKernel(m_context->command_queue(), m_kernel_Daxpy, 1, nullptr, &num,
             nullptr, 0, nullptr, &run_event);

    clWaitForEvents(1, &run_event);
}
