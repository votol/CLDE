#pragma once
#include "IDEOperator.h"
#include <list>
#include <vector>
#include <map>
#include <optional>
#include <complex>

namespace clde
{
struct Monomial
{
    unsigned int outInd = 0;
    double coe = 0.0;
    std::list<unsigned int> inInds;
    std::optional<unsigned int> tFunc;
};

struct MonomialC
{
    unsigned int outInd = 0;
    std::complex<double> coe;
    std::list<unsigned int> inInds;
    std::list<unsigned int> inIndsC;
    std::optional<unsigned int> tFunc;
};

std::list<Monomial> convertMonomials(const std::list<MonomialC>&);


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

        std::vector<std::map<std::list<unsigned int>, double, Key>> m_elements;
        std::vector<std::map<std::pair<std::list<unsigned int>, unsigned int>
                        , double, KeyNoise>> m_elements_noise;

        size_t m_size;
        unsigned int m_num;
        unsigned int m_noise_num;
        std::optional<std::list<unsigned int> > m_dims;

        void print(void);
    public:
        DataAligner(const size_t& dim, const unsigned int& num);
        void add(const Monomial&);
        std::list<unsigned int> getDims(void);
        std::vector<double> getCoes(void);
        std::vector<unsigned int> getIndexes(void);
    };
    size_t m_size;
    unsigned int m_num;
    unsigned int align;
    std::shared_ptr<ICLmanager> m_context;
    cl_kernel m_kernel_oper;
    CLDataStorage<double> m_devie_coes;
    CLDataStorage<unsigned int> m_device_inds;

    void buildKernels(DataAligner& in);
public:
    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    static size_t calculateDimension(_InputIterator begin,  _InputIterator end)
    {
        size_t result = 0;
        for(auto it = begin; it != end; ++it)
        {
            if(it->outInd + 1 > result)
                result = it->outInd + 1;
            for(auto inIt = it->inInds.begin(); inIt != it->inInds.end(); ++inIt)
            {
                if(*inIt + 1 > result)
                    result = *inIt + 1;
            }
        }
        return result;
    }

    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    PolynomialOperator(_InputIterator begin,  _InputIterator end,
                       const unsigned int& num, const std::shared_ptr<ICLmanager>& context ):m_num(num),m_context(context),
                       m_devie_coes(context), m_device_inds(context)
    {
        size_t dim = calculateDimension(begin,end);
        m_size = dim * num + 1;
        DataAligner aligner(dim, num);
        for(auto it = begin; it != end; ++it)
        {
            aligner.add(*it);
        }
        buildKernels(aligner);
    }
    const size_t & dimension() override{return m_size;}
    void apply(const CLDataStorage<double> &in, CLDataStorage<double> &out, const std::vector<double> &) override;
};
};
