#pragma once
#include "IDEOperator.h"
#include "FakeFuncCalculator.h"
#include "PolynomialUtils.h"
#include <list>
#include <vector>
#include <map>
#include <optional>
#include <complex>
#include <iostream>

namespace clde
{
class PolynomialOperator: public IDEOperator
{
private:
    class DataAligner
    {
    private:
        struct Key
        {
            bool operator()(const std::list<unsigned int>&, const std::list<unsigned int>&) const;
        };
        struct KeyNoise
        {
            bool operator()(const std::pair<std::list<unsigned int>, unsigned int>&,
                            const std::pair<std::list<unsigned int>, unsigned int>&) const;
        };

        struct OutElementData
        {
            std::map<std::list<unsigned int>, double, Key> constElements;
            std::map<std::pair<std::list<unsigned int>, unsigned int>
                            , double, KeyNoise> timeElements;
        };

        std::vector<OutElementData> m_elements;

        OperatorDimension m_size;
        unsigned int m_num;
        unsigned int m_time_num;
        std::optional<std::list<unsigned int> > m_dims_const;
        std::optional<std::list<unsigned int> > m_dims_time;

    public:
        DataAligner(const OperatorDimension& dim, const unsigned int& num);
        void add(const Monomial&);
        std::list<unsigned int> getConstDims(void);
        std::list<unsigned int> getTimeDims(void);
        std::vector<double> getCoes(void);
        std::vector<unsigned int> getIndexes(void);
        std::vector<unsigned int> getTimeIndexes(void);
        unsigned int getTimeFuncNumber(void);
    };
    OperatorDimension m_size;
    unsigned int m_num;
    unsigned int m_time_num;
    bool m_is_constant = true;
    std::shared_ptr<ICLmanager> m_context;
    cl_kernel m_kernel_oper;
    size_t m_global_work_size;
    size_t m_local_work_size;
    CLDataStorage<double> m_devie_coes;
    CLDataStorage<unsigned int> m_device_inds;
    CLDataStorage<unsigned int> m_device_time_inds;
    std::shared_ptr<IFuncCalculator> m_time_calc;



    void buildKernels(DataAligner& in);
public:
    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    static OperatorDimension calculateDimension(_InputIterator begin,  _InputIterator end, const bool& square)
    {
        OperatorDimension result;
        for(auto it = begin; it != end; ++it)
        {
            if(it->outInd + 1 > result.out_dim)
                result.out_dim = it->outInd + 1;
            for(auto inIt = it->inInds.begin(); inIt != it->inInds.end(); ++inIt)
            {
                if(*inIt + 1 > result.in_dim)
                    result.in_dim = *inIt + 1;
            }
        }
        if(square)
        {
            size_t size_tmp = (result.in_dim > result.out_dim)? result.in_dim: result.out_dim;
            result.in_dim = size_tmp;
            result.out_dim = size_tmp;
        }
        return result;
    }

    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    PolynomialOperator(_InputIterator begin,  _InputIterator end,
                       const unsigned int& num, const OperatorDimension& dim,
                       const std::shared_ptr<ICLmanager>& context):m_num(num),m_context(context),
                       m_devie_coes(context), m_device_inds(context), m_device_time_inds(context)
    {
        m_size.in_dim = dim.in_dim * num + 1;
        m_size.out_dim = dim.out_dim * num + 1;
        DataAligner aligner(dim, num);
        for(auto it = begin; it != end; ++it)
        {
            aligner.add(*it);
        }
        buildKernels(aligner);

        m_local_work_size = 0;
        cl_int err = 0;
        err = clGetKernelWorkGroupInfo(m_kernel_oper,
                                       m_context->device(),
                                       CL_KERNEL_WORK_GROUP_SIZE,
                                       sizeof(size_t),
                                       &m_local_work_size,
                                       nullptr);

        if(err < 0) {
            throw std::runtime_error("PolynomialOperator_build : OpenCL: problems with kernels work group size");
        }
        m_local_work_size = std::min(m_local_work_size, static_cast<size_t>(64));
        if(m_size.out_dim % m_local_work_size == 0)
            m_global_work_size = m_size.out_dim;
        else
            m_global_work_size = (m_size.out_dim / m_local_work_size + 1) * m_local_work_size;
    }
    ~PolynomialOperator() override;
    const OperatorDimension & dimension() override{return m_size;}
    void apply(const CLDataStorage<double> &in, CLDataStorage<double> &out, const std::vector<double> &) override;
    void setTimeFuncCalculator(const std::shared_ptr<IFuncCalculator>& calc);
};
};
