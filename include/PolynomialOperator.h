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
        struct Comparator
        {
            bool operator() (const std::list<unsigned int>& lhs, const std::list<unsigned int>& rhs);
        };


    public:
        DataAligner(const size_t& dim);
        void add(const Monomial&);
        const size_t& size() const;

    };
    size_t m_size;
    std::shared_ptr<ICLmanager> m_context;

    void buildKernels(const DataAligner& in);
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
    PolynomialOperator(_InputIterator begin,  _InputIterator end, const unsigned int& )
    {
        DataAligner aligner(calculateDimension(begin,end));
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
