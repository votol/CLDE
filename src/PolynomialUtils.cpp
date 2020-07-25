#include "PolynomialUtils.h"
#include <memory>
#include <set>
#include <map>
#include <vector>
#include <iostream>

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
            return false;
        else if(!in.tFunc)
            return true;
        else
            return *tFunc < *in.tFunc;
    }
    void print() const
    {
        std::cout << *coe << " : ";
        std::cout << outInd << " : ";
        for(auto it = inInds.begin(); it != inInds.end(); ++it)
            std::cout << *it << ", ";
        std::cout << ";";
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

std::ostream& clde::operator<<(std::ostream& os, const Monomial& mo)
{
    os << "coefficient: "<< mo.coe;
    os << " output: " <<mo.outInd;
    os << " inds: ";
    for(auto it = mo.inInds.begin(); it != mo.inInds.end(); ++it)
    {
        os << *it;
        if( it != --mo.inInds.end())
            os << ", ";
    }
    os << "; ";
    return os;
}

std::ostream& clde::operator<<(std::ostream& os, const MonomialC& mo)
{
    os << "coefficient: ";
    if(mo.coe.real() != 0.0)
    {
        os << mo.coe.real();
        if(mo.coe.imag() != 0.0)
        {
            if(mo.coe.imag() > 0.0)
                os << " + ";
            os << mo.coe.imag();
            os << "i";
        }
    }
    else
    {
        if(mo.coe.imag() != 0.0)
        {
            os << mo.coe.imag();
            os << "i";
        }
        else
            os << "0.0";
    }
    os << " output: " <<mo.outInd;
    os << " inds: ";
    for(auto it = mo.inInds.begin(); it != mo.inInds.end(); ++it)
    {
        os << *it;
        if( it != --mo.inInds.end())
            os << ", ";
    }
    os << " indsc: ";
    for(auto it = mo.inIndsC.begin(); it != mo.inIndsC.end(); ++it)
    {
        os << *it;
        if( it != --mo.inIndsC.end())
            os << ", ";
    }

    os << "; ";
    return os;
}

std::ofstream& clde::operator<<(std::ofstream& os, const Polynomial& poly)
{
    unsigned int out_int;
    double out_double;

    //number of monomials
    out_int = static_cast<unsigned int>(poly.size());
    os.write(reinterpret_cast <char*>(&out_int), 4);
    for(auto&& mon: poly)
    {
        out_double = mon.coe;
        os.write(reinterpret_cast <char*>(&out_double), 8);
        out_int = mon.outInd;
        os.write(reinterpret_cast <char*>(&out_int), 4);

        out_int = static_cast<unsigned int>(mon.inInds.size());
        os.write(reinterpret_cast <char*>(&out_int), 4);

        for(auto&& ind: mon.inInds)
        {
            out_int = ind;
            os.write(reinterpret_cast <char*>(&out_int), 4);
        }
    }

    return os;
}

std::ifstream& clde::operator>>(std::ifstream& is, Polynomial& poly)
{
    unsigned int mon_count;
    unsigned int inds_count;
    unsigned int ind_read;
    double out_double;

    //number of monomials
    is.read(reinterpret_cast <char*>(&mon_count), 4);
    for(unsigned int mon_ind = 0; mon_ind < mon_count; ++mon_ind)
    {
        poly.push_back(Monomial());
        is.read(reinterpret_cast <char*>(&out_double), 8);
        poly.back().coe = out_double;

        is.read(reinterpret_cast <char*>(&ind_read), 4);
        poly.back().outInd = ind_read;

        is.read(reinterpret_cast <char*>(&inds_count), 4);
        for(unsigned int ind = 0; ind < inds_count; ++ind)
        {
            is.read(reinterpret_cast <char*>(&ind_read), 4);
            poly.back().inInds.push_back(ind_read);
        }
    }

    return is;
}

void clde::unitePoly(PolynomialC& dest, PolynomialC&& source, const std::complex<double>& coe)
{
    for(auto&& mon: source)
    {
        mon.coe *= coe;
    }
    dest.splice(dest.end(), source);
}

void clde::unitePoly(Polynomial& dest, Polynomial&& source, const double& coe)
{
    for(auto&& mon: source)
    {
        mon.coe *= coe;
    }
    dest.splice(dest.end(), source);
}

Polynomial clde::convertMonomials(const PolynomialC& in)
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

    Polynomial result;
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

Polynomial clde::polynomialDerivative(const Polynomial& in, const unsigned int& var)
{
    std::set<MonomialAnc> monomials;
    for(auto in_it = in.begin(); in_it != in.end(); ++ in_it)
    {
        auto inds_tmp = in_it->inInds;
        inds_tmp.sort();
        unsigned int count = 0;
        for(auto ind_it = inds_tmp.begin(); ind_it != inds_tmp.end(); ++ind_it)
        {
            if(*ind_it == var)
            {
                MonomialAnc element_ins;
                element_ins.outInd = in_it->outInd;
                *(element_ins.coe) = in_it->coe;
                element_ins.inInds = inds_tmp;
                auto erase_it = element_ins.inInds.begin();
                std::advance(erase_it, count);
                element_ins.inInds.erase(erase_it);
                element_ins.tFunc = in_it->tFunc;
                auto adding_it = monomials.find(element_ins);
                if(adding_it != monomials.end())
                    *(adding_it->coe) += *(element_ins.coe);
                else {
                    monomials.insert(element_ins);
                }
            }
            count++;
        }
    }
    Polynomial result;
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

std::vector<Polynomial> clde::polynomialGradient(const Polynomial& in, const unsigned int& var_count)
{
    std::vector<Polynomial> result(var_count);
    for(unsigned int ind = 0; ind < var_count; ++ind)
        result[ind] = polynomialDerivative(in, ind);
    return result;
}
