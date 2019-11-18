#include "PolynomialOperator.h"
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <map>
#include <vector>
#include "CLBLAS.h"

using namespace clde;

bool PolynomialOperator::DataAligner::Key::operator()(const std::list<unsigned int> &lhs, const std::list<unsigned int> &rhs) const
{
    if(lhs.size() != rhs.size())
        return lhs.size()<rhs.size();
    auto lhs_it = lhs.begin();
    auto rhs_it = rhs.begin();
    for(;lhs_it != lhs.end() && rhs_it != rhs.end(); ++lhs_it, ++rhs_it)
    {
        if(*lhs_it != *rhs_it)
            return *lhs_it < *rhs_it;
    }
    return false;
}

bool PolynomialOperator::DataAligner::KeyNoise::operator()(const std::pair<std::list<unsigned int>,
                                                           unsigned int> &lhs, const std::pair<std::list<unsigned int>, unsigned int> &rhs)const
{
    if(lhs.first.size() != rhs.first.size())
        return lhs.first.size()<rhs.first.size();
    auto lhs_it = lhs.first.begin();
    auto rhs_it = rhs.first.begin();
    for(;lhs_it != lhs.first.end() && rhs_it != rhs.first.end(); ++lhs_it, ++rhs_it)
    {
        if(*lhs_it != *rhs_it)
            return *lhs_it < *rhs_it;
    }
    return lhs.second < rhs.second;
}

PolynomialOperator::DataAligner::DataAligner(const OperatorDimension& dim, const unsigned int& num):
    m_elements(dim.out_dim + 1), m_size(dim), m_num(num), m_time_num(0)
{

}

void PolynomialOperator::DataAligner::add(const Monomial &in)
{
    std::list<unsigned int> inInds = in.inInds;
    inInds.sort();
    for(auto it = inInds.begin(); it != inInds.end(); ++it)
    {
        *it += 1;
    }
    if(in.tFunc)
    {
        if(*(in.tFunc) + 1 > m_time_num)
            m_time_num = *(in.tFunc) + 1;
        std::pair<std::list<unsigned int>, unsigned int> ind(inInds, *(in.tFunc));
        std::map<std::pair<std::list<unsigned int>, unsigned int>
              , double, KeyNoise>& output = m_elements[in.outInd + 1].timeElements;
        auto it = output.find(ind);
        if(it != output.end())
        {
            it->second += in.coe;
        }
        else
        {
            output.insert(std::pair<std::pair<std::list<unsigned int>, unsigned int>
                          , double>(ind, in.coe));
        }
    }
    else
    {
        std::map<std::list<unsigned int>, double, Key>& output = m_elements[in.outInd + 1].constElements;
        auto it = output.find(inInds);
        if(it != output.end())
        {
            it->second += in.coe;
        }
        else
        {
            output.insert(std::pair<std::list<unsigned int>, double>(inInds, in.coe));
        }
    }
    if(m_dims_const)
        m_dims_const.reset();
    if(m_dims_time)
        m_dims_time.reset();
}

std::list<unsigned int> PolynomialOperator::DataAligner::getConstDims()
{
    if(m_dims_const)
        return *m_dims_const;
    std::list<unsigned int> result;

    for(auto outIt = m_elements.begin();outIt != m_elements.end();++outIt)
    {
        int adding_size = (static_cast<int>(outIt->constElements.size()) - static_cast<int>(result.size()));
        for(int ind = 0; ind< adding_size; ++ind)
        {
            result.push_front(0);
        }
        auto resIt = result.rbegin();
        auto inIt = outIt->constElements.rbegin();
        for(; inIt != outIt->constElements.rend();++inIt, ++resIt)
        {
            if(inIt->first.size() > *resIt)
                *resIt = static_cast<unsigned int>(inIt->first.size());
        }
    }

    m_dims_const = result;
    return result;
}

std::list<unsigned int> PolynomialOperator::DataAligner::getTimeDims()
{
    if(m_dims_time)
        return *m_dims_time;
    std::list<unsigned int> result;

    for(auto outIt = m_elements.begin();outIt != m_elements.end();++outIt)
    {
        int adding_size = (static_cast<int>(outIt->timeElements.size()) - static_cast<int>(result.size()));
        for(int ind = 0; ind< adding_size; ++ind)
        {
            result.push_front(0);
        }
        auto resIt = result.rbegin();
        auto inIt = outIt->timeElements.rbegin();
        for(; inIt != outIt->timeElements.rend();++inIt, ++resIt)
        {
            if(inIt->first.first.size() > *resIt)
                *resIt = static_cast<unsigned int>(inIt->first.first.size());
        }
    }

    m_dims_time = result;
    return result;
}

std::vector<double> PolynomialOperator::DataAligner::getCoes()
{
    if(!m_dims_const)
        getConstDims();
    if(!m_dims_time)
        getTimeDims();

    std::vector<double> result((m_dims_const->size() + m_dims_time->size())* (1 + m_num * m_size.out_dim));
    auto resIt = result.begin();
    for(auto dimIt = m_dims_const->begin();dimIt != m_dims_const->end();++dimIt)
    {
        *resIt = 0.0;
        ++resIt;
    }
    for(auto dimIt = m_dims_time->begin();dimIt != m_dims_time->end();++dimIt)
    {
        *resIt = 0.0;
        ++resIt;
    }
    for(unsigned int ind = 0;ind<m_num;++ind)
    {
        for(auto outIt = ++(m_elements.begin());outIt != m_elements.end();++outIt)
        {
            auto constIt = outIt->constElements.rbegin();
            for(auto dimIt = m_dims_const->rbegin();dimIt != m_dims_const->rend();++dimIt)
            {
                if(constIt == outIt->constElements.rend())
                    *resIt = 0.0;
                else
                {
                    *resIt = constIt->second;
                    ++constIt;
                }
                ++resIt;
            }

            auto timeIt = outIt->timeElements.rbegin();
            for(auto dimIt = m_dims_time->rbegin();dimIt != m_dims_time->rend();++dimIt)
            {
                if(timeIt == outIt->timeElements.rend())
                    *resIt = 0.0;
                else
                {
                    *resIt = timeIt->second;
                    ++timeIt;
                }
                ++resIt;
            }
        }
    }
    return result;
}

std::vector<unsigned int> PolynomialOperator::DataAligner::getIndexes()
{
    if(!m_dims_const)
        getConstDims();
    if(!m_dims_time)
        getTimeDims();
    unsigned int one_out_size_const = 0;
    for(auto it = m_dims_const->begin();it != m_dims_const->end(); ++it)
    {
        one_out_size_const += *it;
    }
    unsigned int one_out_size_time = 0;
    for(auto it = m_dims_time->begin();it != m_dims_time->end(); ++it)
    {
        one_out_size_time += *it;
    }

    std::vector<unsigned int> result((one_out_size_const + one_out_size_time) * (1+ m_num * m_size.out_dim));
    auto resIt = result.begin();
    for(unsigned int ind = 0; ind < one_out_size_const + one_out_size_time;++ind)
    {
        *resIt = 0;
        ++resIt;
    }

    for(unsigned int ind = 0; ind<m_num; ++ind)
    {
        for(auto outIt = ++(m_elements.begin());outIt != m_elements.end();++outIt)
        {
            auto control_indexs_number_it = m_dims_const->rbegin();
            unsigned int index_added = 0;
            for(auto summandIt = outIt->constElements.rbegin();summandIt != outIt->constElements.rend();++summandIt)
            {
                unsigned int single_index_added = 0;
                for(auto indexIt = summandIt->first.begin(); indexIt != summandIt->first.end(); indexIt++)
                {
                    *resIt = *indexIt + ind * static_cast<unsigned int>(m_size.in_dim);
                    ++resIt;
                    ++single_index_added;
                }
                for(;single_index_added<*control_indexs_number_it;++single_index_added)
                {
                    *resIt = 0;
                    ++resIt;
                }
                ++control_indexs_number_it;
                index_added += single_index_added;
            }
            for(;index_added<one_out_size_const;++index_added)
            {
                *resIt = 0;
                ++resIt;
            }

            control_indexs_number_it = m_dims_time->rbegin();
            index_added = 0;
            for(auto summandIt = outIt->timeElements.rbegin();summandIt != outIt->timeElements.rend();++summandIt)
            {
                unsigned int single_index_added = 0;
                for(auto indexIt = summandIt->first.first.begin(); indexIt != summandIt->first.first.end(); indexIt++)
                {
                    *resIt = *indexIt + ind * static_cast<unsigned int>(m_size.in_dim);
                    ++resIt;
                    ++single_index_added;
                }
                for(;single_index_added<*control_indexs_number_it;++single_index_added)
                {
                    *resIt = 0;
                    ++resIt;
                }
                ++control_indexs_number_it;
                index_added += single_index_added;

            }
            for(;index_added<one_out_size_time;++index_added)
            {
                *resIt = 0;
                ++resIt;
            }
        }
    }
    return result;
}

std::vector<unsigned int> PolynomialOperator::DataAligner::getTimeIndexes()
{
    if(!m_dims_time)
        getTimeDims();
    std::vector<unsigned int> result(m_dims_time->size() * (1+ m_num * m_size.out_dim));
    auto resIt = result.begin();
    for(auto dimIt = m_dims_time->begin();dimIt != m_dims_time->end();++dimIt)
    {
        *resIt = 0;
        ++resIt;
    }
    for(unsigned int ind = 0; ind<m_num; ++ind)
    {
        for(auto outIt = ++(m_elements.begin());outIt != m_elements.end();++outIt)
        {
            auto timeIt = outIt->timeElements.rbegin();
            for(auto dimIt = m_dims_time->rbegin();dimIt != m_dims_time->rend();++dimIt)
            {
                if(timeIt == outIt->timeElements.rend())
                    *resIt = 0;
                else
                {
                    *resIt = timeIt->first.second + ind * m_time_num;
                    ++timeIt;
                }
                ++resIt;
            }
        }
    }
    return result;
}

unsigned int PolynomialOperator::DataAligner::getTimeFuncNumber()
{
    return m_time_num * m_num;
}

void PolynomialOperator::buildKernels(DataAligner &in)
{
    std::list<unsigned int> dims_const(in.getConstDims());
    std::list<unsigned int> dims_time(in.getTimeDims());
    unsigned int one_out_size = 0;
    for(auto it = dims_const.begin();it != dims_const.end(); ++it)
    {
        one_out_size += *it;
    }
    for(auto it = dims_time.begin();it != dims_time.end(); ++it)
    {
        one_out_size += *it;
    }

    /*std::cout<<std::endl<<"--------------------------"<<std::endl;
    std::vector<double> coes(in.getCoes());
    std::vector<unsigned int> inInds(in.getIndexes());
    auto coesIt = coes.begin();
    auto inIndsIt = inInds.begin();
    for(unsigned int ind = 0; ind<m_size.out_dim;++ind)
    {
        for(auto it = dims_const.rbegin();it != dims_const.rend(); ++it)
        {
            std::cout << *coesIt << " , ";
            ++coesIt;
            for (unsigned ind = 0; ind < *it; ++ind)
            {
                std::cout << *inIndsIt;
                if(ind != *it -1)
                    std::cout << " , ";
                ++inIndsIt;
            }
            std::cout << " ; ";
        }
        for(auto it = dims_time.rbegin();it != dims_time.rend(); ++it)
        {
            std::cout << *coesIt << " , ";
            ++coesIt;
            for (unsigned ind = 0; ind < *it; ++ind)
            {
                std::cout << *inIndsIt;
                if(ind != *it -1)
                    std::cout << " , ";
                ++inIndsIt;
            }
            std::cout << " ; ";
        }
        std::cout << std::endl;
    }*/

    m_devie_coes = in.getCoes();
    m_device_inds = in.getIndexes();
    if(dims_time.size() != 0)
    {
        m_is_constant = false;
        m_device_time_inds = in.getTimeIndexes();
        m_time_calc = std::shared_ptr<IFuncCalculator>(new FakeFuncCalculator(m_context));
        m_time_calc->init(in.getTimeFuncNumber());
        m_time_num = in.getTimeFuncNumber();
    }


    std::string kernel_str;
    kernel_str += "__kernel void Oper(__global double* out , __global double* in, __global double* coes, __global unsigned int* inds";
    if(dims_time.size() == 0)
    {
        kernel_str += ") {\n";
    }
    else {
        kernel_str += ", __global double* time_coes, const unsigned int time_offset, __global unsigned int* time_inds) {\n";
    }
    kernel_str += "unsigned int gid = get_global_id(0);\n";
    kernel_str += "__global double* loc_coes = coes ;\n";
    kernel_str += "loc_coes += gid * ";
    kernel_str += std::to_string(dims_const.size() + dims_time.size());
    kernel_str += ";\n";
    kernel_str += "__global unsigned int* loc_inds = inds ;\n";
    kernel_str += "loc_inds += gid * ";
    kernel_str += std::to_string(one_out_size);
    kernel_str += ";\n";
    kernel_str += "out[gid] = 0.0;\n";
    kernel_str += "double result;\n";
    for(auto it = dims_const.rbegin(); it != dims_const.rend(); ++it)
    {
        kernel_str += "result = *loc_coes;\n";
        kernel_str += "++loc_coes;\n";
        for(unsigned int ind = 0; ind < *it; ++ind)
        {
            kernel_str += "result *= in[*loc_inds];\n";
            kernel_str += "++loc_inds;\n";
        }
        kernel_str += "out[gid] += result;\n";
    }
    if(dims_time.size() != 0)
    {
        kernel_str += "__global double* loc_time_coes = time_coes ;\n";
        kernel_str += "loc_time_coes += time_offset;\n";
        kernel_str += "__global unsigned int* loc_time_inds = time_inds ;\n";
        kernel_str += "loc_time_inds += gid * ";
        kernel_str += std::to_string(dims_time.size());
        kernel_str += ";\n";
        for(auto it = dims_time.rbegin(); it != dims_time.rend(); ++it)
        {
            kernel_str += "result = *loc_coes;\n";
            kernel_str += "++loc_coes;\n";
            kernel_str += "result *= loc_time_coes[*loc_time_inds];\n";
            kernel_str += "++loc_time_inds;\n";
            for(unsigned int ind = 0; ind < *it; ++ind)
            {
                kernel_str += "result *= in[*loc_inds];\n";
                kernel_str += "++loc_inds;\n";
            }
            kernel_str += "out[gid] += result;\n";
        }
    }
    kernel_str += "}";
    //std::cout << kernel_str << std::endl;
    m_kernel_oper = CLBLAS::createKernel(kernel_str.c_str(), "Oper", m_context);
    cl_int err = 0;
    err = clSetKernelArg(m_kernel_oper, 2, sizeof(cl_mem), &m_devie_coes.data());
    err |= clSetKernelArg(m_kernel_oper, 3, sizeof(cl_mem), &m_device_inds.data());
    if(dims_time.size() != 0)
    {
        err |= clSetKernelArg(m_kernel_oper, 4, sizeof(cl_mem), &m_time_calc->data().data());
        err |= clSetKernelArg(m_kernel_oper, 6, sizeof(cl_mem), &m_device_time_inds.data());
    }
    if(err < 0) {
        throw std::runtime_error("PolynomialOperator_init : OpenCL: couldn't set an argument for the polynomial operator kernel");
    }
}

PolynomialOperator::~PolynomialOperator()
{
    clReleaseKernel(m_kernel_oper);
}

void PolynomialOperator::apply(const CLDataStorage<double> &in, CLDataStorage<double> &out, const std::vector<double> &param)
{
    if(in.size() != m_size.in_dim && out.size() != m_size.out_dim)
        throw std::runtime_error("Polynomial Operator: vectors size does not agree with operator dimensions");

    size_t num = out.size();
    cl_int err = 0;
    cl_event run_event;

    err = clSetKernelArg(m_kernel_oper, 0, sizeof(cl_mem), &out.data());
    err |= clSetKernelArg(m_kernel_oper, 1, sizeof(cl_mem), &in.data());
    if(!m_is_constant)
    {
        err |= clSetKernelArg(m_kernel_oper, 5, sizeof(unsigned int), &m_time_calc->process(param[0]));
    }

    if(err < 0) {
        throw std::runtime_error("PolynomialOperator_run : OpenCL: couldn't set an argument for the operator kernel");
    }

    err = clEnqueueNDRangeKernel(m_context->command_queue(), m_kernel_oper, 1, nullptr, &num,
             nullptr, 0, nullptr, &run_event);

    if(err < 0) {
        throw std::runtime_error("PolynomialOperator_run : OpenCL: problems with kernels");
    }
    clWaitForEvents(1, &run_event);
    /*cl_int status;
    err = clGetEventInfo(run_event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &status, nullptr);
    if(err < 0 )
    {
        throw std::runtime_error("PolynomialOperator_run : OpenCL: can't get status");
    }
    if(status < 0)
    {
        throw std::runtime_error("PolynomialOperator_run : OpenCL: kernel faild");
    }*/
    clReleaseEvent(run_event);
}

void PolynomialOperator::setTimeFuncCalculator(const std::shared_ptr<IFuncCalculator> &calc)
{
    if(m_is_constant)
        return;
    cl_int err = 0;
    m_time_calc = calc;
    m_time_calc->init(m_time_num);
    err = clSetKernelArg(m_kernel_oper, 4, sizeof(cl_mem), &m_time_calc->data().data());
    if(err < 0) {
        throw std::runtime_error("OpenCL: couldn't set an argument for the polynomial operator kernel");
    }
}
