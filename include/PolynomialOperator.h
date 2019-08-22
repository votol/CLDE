#pragma once
#include "IDEOperator.h"
#include <list>
#include <vector>
#include <map>
#include <optional>

namespace clde
{
struct Monomial
{
    unsigned int outInd = 0;
    double coe = 0.0;
    std::list<unsigned int> inInds;
    std::optional<unsigned int> noise;
};


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

        bool m_aligned;
        size_t m_size;
        unsigned int m_num;
        unsigned int m_noise_num;

        void print(void);
    public:
        DataAligner(const size_t& dim, const unsigned int& num);
        void add(const Monomial&);
        void align(void);
    };
    size_t m_size;
    std::shared_ptr<ICLmanager> m_context;

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
                       const unsigned int& num, const std::shared_ptr<ICLmanager>& context ):m_context(context)
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
    void apply(const CLDataStorage<double> &in, CLDataStorage<double> &out, const std::vector<double> &) override
    {
        out = in;
    }
};
};
