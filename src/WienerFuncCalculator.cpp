#include "WienerFuncCalculator.h"
#include "kernels/WienerCalc.h"
#include "CLBLAS.h"
#include <string>
#include <random>
#include <stdexcept>

using namespace clde;

WienerFuncCalculator::WienerFuncCalculator(const std::shared_ptr<ICLmanager>& context):
    m_context(context), m_data(context), m_random_numbers(context)
{

}
WienerFuncCalculator::~WienerFuncCalculator()
{
    if(m_kernel)
        clReleaseKernel(*m_kernel);
}

void WienerFuncCalculator::init(const size_t& num)
{
    m_threads_num = 32 * (num/32);
    m_threads_num += (m_threads_num != num)? 32: 0;
    std::string kernel;
    kernel = "#define DIM " + std::to_string(num) + "\n";
    kernel += WienerCalc_kernel;
    if(m_kernel)
        clReleaseKernel(*m_kernel);
    m_kernel = CLBLAS::createKernel(kernel.c_str(), "wienerCalc", m_context);
    m_data = std::vector<double>(num, 0.0);

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<double> init_rand(97 * (m_threads_num/32));
    for(auto it = init_rand.begin(); it != init_rand.end(); ++it)
    {
        *it = dis(gen);
    }
    m_random_numbers = init_rand;
    m_random_shift = 0;

    cl_int err = 0;
    err = clSetKernelArg(*m_kernel, 0, sizeof(cl_mem), &m_data.data());
    err |= clSetKernelArg(*m_kernel, 1, sizeof(cl_mem), &m_random_numbers.data());
    if(err < 0) {
        throw std::runtime_error("Wiener calculator init : OpenCL: couldn't set an argument for the operator kernel");
    }

}
const CLDataStorage<double>& WienerFuncCalculator::data() const
{
    return m_data;
}
const unsigned int& WienerFuncCalculator::process(const double& in_t)
{
    if(in_t <= m_current_time)
        return m_ind;
    cl_event run_event;
    cl_int err = 0;
    m_delta_time = in_t - m_current_time;
    err = clSetKernelArg(*m_kernel, 2, sizeof(unsigned int), &m_random_shift);
    err = clSetKernelArg(*m_kernel, 3, sizeof(double), &m_delta_time);
    if(err < 0) {
        throw std::runtime_error("Wiener calculator process : OpenCL: couldn't set an argument for the operator kernel");
    }

    err = clEnqueueNDRangeKernel(m_context->command_queue(), *m_kernel, 1, nullptr, &m_threads_num,
             &m_work_group, 0, nullptr, &run_event);

    clWaitForEvents(1, &run_event);
    clReleaseEvent(run_event);

    m_random_shift = (m_random_shift + 32) % 97;
    m_current_time = in_t;
    return m_ind;
}
