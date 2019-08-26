#include "PolynomialOperator.h"
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <map>
#include <vector>
#include "CLBLAS.h"

using namespace clde;

struct MonomialAnc
{
    unsigned int outInd = 0;
    std::shared_ptr<double> coe;
    std::list<unsigned int> inInds;
    std::optional<unsigned int> tFunc;
    MonomialAnc():coe(new double){*coe = 0.0;}
    bool operator<(const MonomialAnc& in) const
    {
        if(outInd != in.outInd)
            return outInd < in.outInd;
        auto m_it = inInds.begin();
        auto in_it = in.inInds.begin();
        for(;m_it != inInds.end() && in_it != in.inInds.end(); ++ m_it, ++ in_it)
        {
            if(*m_it != *in_it)
                return *m_it < *in_it;
        }
        if(inInds.size() != in.inInds.size())
            return inInds.size() < in.inInds.size();
        if(!tFunc)
            return true;
        else if(!in.tFunc)
            return false;
        else
            return *tFunc < *in.tFunc;
    }
};

class MonomialGenerator
{
private:
    std::list<unsigned int> m_norm;
    std::list<unsigned int> m_conj;
    std::vector<bool> m_state;
    bool m_overloaded = false;
public:
    MonomialGenerator(const std::list<unsigned int>& norm, const std::list<unsigned int>& conj):
        m_norm(norm), m_conj(conj), m_state(norm.size() + conj.size(), false)
    {

    }
    MonomialGenerator& operator++()
    {
        if(m_overloaded)
            return *this;

        auto it = m_state.begin();
        for(;;)
        {
            if(it == m_state.end())
            {
                m_overloaded = true;
                break;
            }
            if(*it)
            {
                *it = false;
            }
            else
            {
                *it = true;
                break;
            }

            ++it;
        }
        return *this;
    }
    bool finished(void)
    {
        return m_overloaded;
    }
    std::pair<std::list<unsigned int>, int> generate(void)
    {
        std::list<unsigned int> result;
        int coe = 0;
        if(m_overloaded)
            return std::pair<std::list<unsigned int>, int>(result,coe);

        auto it = m_state.begin();
        for(auto it_norm = m_norm.begin();it_norm != m_norm.end();++it_norm, ++it)
        {
            if(*it)
            {
                coe++;
                result.push_back(2 * (*it_norm) + 1);
            }
            else
            {
                result.push_back(2 * (*it_norm));
            }
        }
        for(auto it_conj = m_conj.begin();it_conj != m_conj.end();++it_conj, ++it)
        {
            if(*it)
            {
                coe--;
                result.push_back(2 * (*it_conj) + 1);
            }
            else
            {
                result.push_back(2 * (*it_conj));
            }
        }
        return std::pair<std::list<unsigned int>, int>(result,coe);
    }
};

struct KeyListComp
{
    bool operator()(const std::list<unsigned int> &lhs, const std::list<unsigned int> &rhs) const
    {
        auto lhs_it = lhs.begin();
        auto rhs_it = rhs.begin();
        for(;lhs_it != lhs.end() && rhs_it != rhs.end(); ++ lhs_it, ++ rhs_it)
        {
            if(*lhs_it != *rhs_it)
                return *lhs_it < *rhs_it;
        }
        return lhs.size() < rhs.size();
    }
};

std::list<Monomial> clde::convertMonomials(const std::list<MonomialC>& in)
{
    std::set<MonomialAnc> monomials;
    for(auto in_it = in.begin(); in_it != in.end(); ++ in_it)
    {
        std::map<std::list<unsigned int>, int, KeyListComp> mon_real;
        std::map<std::list<unsigned int>, int, KeyListComp> mon_imag;
        MonomialGenerator generator(in_it->inInds, in_it->inIndsC);
        while(!generator.finished())
        {
            auto entry = generator.generate();
            entry.first.sort();
            std::map<std::list<unsigned int>, int, KeyListComp>::iterator find_it;
            if(entry.second < 0)
                entry.second += 4 * (-entry.second/4 + 1);
            switch(entry.second % 4)
            {
            case 0:
                find_it = mon_real.find(entry.first);
                if(find_it != mon_real.end())
                    find_it->second += 1;
                else
                {
                    entry.second = 1;
                    mon_real.insert(entry);
                }
                break;
            case 1:
                find_it = mon_imag.find(entry.first);
                if(find_it != mon_imag.end())
                    find_it->second += 1;
                else
                {
                    entry.second = 1;
                    mon_imag.insert(entry);
                }
                break;
            case 2:
                find_it = mon_real.find(entry.first);
                if(find_it != mon_real.end())
                    find_it->second -= 1;
                else
                {
                    entry.second = -1;
                    mon_real.insert(entry);
                }
                break;
            case 3:
                find_it = mon_imag.find(entry.first);
                if(find_it != mon_imag.end())
                    find_it->second -= 1;
                else
                {
                    entry.second = -1;
                    mon_imag.insert(entry);
                }
                break;
            }
            ++generator;
        }
        for(auto it = mon_real.begin(); it != mon_real.end(); ++it)
        {
            if(it->second == 0)
                continue;
            MonomialAnc element_r;
            element_r.outInd = in_it->outInd * 2;
            element_r.inInds = it->first;
            *(element_r.coe) = double(it->second) * in_it->coe.real();
            element_r.tFunc = in_it->tFunc;
            auto adding_it = monomials.find(element_r);
            if(adding_it != monomials.end())
                *(adding_it->coe) += *(element_r.coe);
            else {
                monomials.insert(element_r);
            }

            MonomialAnc element_i;
            element_i.outInd = in_it->outInd * 2 + 1;
            element_i.inInds = it->first;
            *(element_i.coe) = double(it->second) * in_it->coe.imag();
            element_i.tFunc = in_it->tFunc;
            adding_it = monomials.find(element_i);
            if(adding_it != monomials.end())
                *(adding_it->coe) += *(element_i.coe);
            else {
                monomials.insert(element_i);
            }
        }

        for(auto it = mon_imag.begin(); it != mon_imag.end(); ++it)
        {
            if(it->second == 0)
                continue;
            MonomialAnc element_r;
            element_r.outInd = in_it->outInd * 2;
            element_r.inInds = it->first;
            *(element_r.coe) = -double(it->second) * in_it->coe.imag();
            element_r.tFunc = in_it->tFunc;
            auto adding_it = monomials.find(element_r);
            if(adding_it != monomials.end())
                *(adding_it->coe) += *(element_r.coe);
            else {
                monomials.insert(element_r);
            }

            MonomialAnc element_i;
            element_i.outInd = in_it->outInd * 2 + 1;
            element_i.inInds = it->first;
            *(element_i.coe) = double(it->second) * in_it->coe.real();
            element_i.tFunc = in_it->tFunc;
            adding_it = monomials.find(element_i);
            if(adding_it != monomials.end())
                *(adding_it->coe) += *(element_i.coe);
            else {
                monomials.insert(element_i);
            }
        }
    }

    std::list<Monomial> result;
    for(auto it = monomials.begin();it != monomials.end(); ++it)
    {
        if(*(it->coe) == 0.0)
            continue;
        Monomial element;
        element.outInd = it->outInd;
        element.inInds = it->inInds;
        element.tFunc = it->tFunc;
        element.coe = *(it->coe);
        result.push_back(element);
    }

    return result;
}

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

PolynomialOperator::DataAligner::DataAligner(const size_t& dim, const unsigned int& num):
    m_elements(dim + 1), m_size(dim), m_num(num), m_time_num(0)
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

    std::vector<double> result((m_dims_const->size() + m_dims_time->size())* (1 + m_num*m_size));
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

    std::vector<unsigned int> result((one_out_size_const + one_out_size_time) * (1+ m_num*m_size));
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
            unsigned int index_added = 0;
            for(auto summandIt = outIt->constElements.rbegin();summandIt != outIt->constElements.rend();++summandIt)
            {
                for(auto indexIt = summandIt->first.begin(); indexIt != summandIt->first.end(); indexIt++)
                {
                    *resIt = *indexIt + ind * static_cast<unsigned int>(m_size);
                    ++resIt;
                    index_added++;
                }

            }
            for(;index_added<one_out_size_const;++index_added)
            {
                *resIt = 0;
                ++resIt;
            }

            index_added = 0;
            for(auto summandIt = outIt->timeElements.rbegin();summandIt != outIt->timeElements.rend();++summandIt)
            {
                for(auto indexIt = summandIt->first.first.begin(); indexIt != summandIt->first.first.end(); indexIt++)
                {
                    *resIt = *indexIt + ind * static_cast<unsigned int>(m_size);
                    ++resIt;
                    index_added++;
                }

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
    std::vector<unsigned int> result(m_dims_time->size() * (1+ m_num*m_size));
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
    for(unsigned int ind = 0; ind<m_size;++ind)
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

void PolynomialOperator::apply(const CLDataStorage<double> &in, CLDataStorage<double> &out, const std::vector<double> &param)
{
    if(in.size() != out.size())
        throw std::runtime_error("Polynomial Operator: input vectors should have the same size");

    size_t num = in.size();
    cl_int err = 0;
    static cl_event run_event;

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

    clWaitForEvents(1, &run_event);
    clReleaseEvent(run_event);
}

void PolynomialOperator::setTimeFuncCalculator(const std::shared_ptr<IFuncCalculator> &calc)
{
    cl_int err = 0;
    m_time_calc = calc;
    m_time_calc->init(m_time_num);
    err = clSetKernelArg(m_kernel_oper, 4, sizeof(cl_mem), &m_time_calc->data());
    if(err < 0) {
        throw std::runtime_error("OpenCL: couldn't set an argument for the polynomial operator kernel");
    }
}
