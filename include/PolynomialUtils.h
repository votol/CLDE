#pragma once
#include <list>
#include <optional>
#include <complex>
#include <vector>
#include <fstream>

namespace clde
{
struct Monomial
{
    unsigned int outInd = 0;
    double coe = 0.0;
    std::list<unsigned int> inInds;
    std::optional<unsigned int> tFunc;
    friend std::ostream& operator<<(std::ostream& os, const Monomial& mo);
};

struct MonomialC
{
    unsigned int outInd = 0;
    std::complex<double> coe;
    std::list<unsigned int> inInds;
    std::list<unsigned int> inIndsC;
    std::optional<unsigned int> tFunc;
};
using Polynomial = std::list<Monomial>;
using PolynomialC = std::list<MonomialC>;


std::ostream& operator<<(std::ostream& os, const Monomial& mo);
std::ofstream& operator<<(std::ofstream& os, const Polynomial& poly);
std::ifstream& operator>>(std::ifstream& os, Polynomial& poly);
Polynomial convertMonomials(const PolynomialC&);
Polynomial polynomialDerivative(const Polynomial&, const unsigned int& var);
std::vector<Polynomial> polynomialGradient(const Polynomial&, const unsigned int& var_count);
}
