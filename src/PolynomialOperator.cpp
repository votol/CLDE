#include "PolynomialOperator.h"
#include <iostream>

using namespace clde;

bool PolynomialOperator::DataAligner::Key::operator()(const std::list<unsigned int> &lhs, const std::list<unsigned int> &rhs) const
{
    auto lhs_it = lhs.begin();
    auto rhs_it = rhs.begin();
    for(;lhs_it != lhs.end() && rhs_it != rhs.end(); ++lhs_it, ++rhs_it)
    {
        if(*lhs_it != *rhs_it)
            return *lhs_it < *rhs_it;
    }
    return lhs.size()<rhs.size();
}

bool PolynomialOperator::DataAligner::KeyNoise::operator()(const std::pair<std::list<unsigned int>,
                                                           unsigned int> &lhs, const std::pair<std::list<unsigned int>, unsigned int> &rhs)const
{
    auto lhs_it = lhs.first.begin();
    auto rhs_it = rhs.first.begin();
    for(;lhs_it != lhs.first.end() && rhs_it != rhs.first.end(); ++lhs_it, ++rhs_it)
    {
        if(*lhs_it != *rhs_it)
            return *lhs_it < *rhs_it;
    }
    if(lhs.first.size() != rhs.first.size())
        return lhs.first.size()<rhs.first.size();
    return lhs.second < rhs.second;
}

PolynomialOperator::DataAligner::DataAligner(const size_t& dim, const unsigned int& num):
    m_elements(dim + 1), m_elements_noise(dim + 1), m_aligned(false), m_size(dim), m_num(num), m_noise_num(0)
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
    if(in.noise)
    {
        if(*(in.noise) + 1 > m_noise_num)
            m_noise_num = *(in.noise) + 1;
        std::pair<std::list<unsigned int>, unsigned int> ind(inInds, *(in.noise));
        std::map<std::pair<std::list<unsigned int>, unsigned int>
              , double, KeyNoise>& output = m_elements_noise[in.outInd + 1];
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
        std::map<std::list<unsigned int>, double, Key>& output = m_elements[in.outInd + 1];
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
    if(m_aligned)
        m_aligned = false;
}

void PolynomialOperator::DataAligner::print()
{
    for(auto it = m_elements.begin(); it != m_elements.end(); ++it)
    {
        for(auto it1 = it->begin();it1!= it->end();++it1)
        {
            std::cout<<it1->second<<",";
            for(auto it2 = it1->first.begin();it2 != it1->first.end(); ++it2)
            {
                std::cout<<*it2<<",";
            }
            std::cout<<";";
        }
        std::cout<<std::endl;
    }
}


void PolynomialOperator::DataAligner::align(void)
{
    if(m_aligned)
        return;

    print();
    m_aligned = true;
}

void PolynomialOperator::buildKernels(DataAligner &in)
{
    in.align();
}
