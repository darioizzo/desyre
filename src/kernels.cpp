// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <dsyre/kernels.hpp>

namespace dsyre
{
// Functions and derivatives
double inv(double a)
{
    return 1. / a;
}
double dinv(double a)
{
    return -1. / a / a;
}
double ddinv(double a)
{
    return 2 / a / a / a;
}
std::string pinv(std::string arg)
{
    return "1/" + arg;
}

double cos(double a)
{
    return std::cos(a);
}
double dcos(double a)
{
    return -std::sin(a);
}
double ddcos(double a)
{
    return -std::cos(a);
}
std::string pcos(std::string arg)
{
    return "cos(" + arg + ")";
}

double sin(double a)
{
    return std::sin(a);
}
double dsin(double a)
{
    return std::cos(a);
}
double ddsin(double a)
{
    return -std::sin(a);
}
std::string psin(std::string arg)
{
    return "sin(" + arg + ")";
}

double exp(double a)
{
    return std::exp(a);
}
double dexp(double a)
{
    return std::exp(a);
}
double ddexp(double a)
{
    return std::exp(a);
}
std::string pexp(std::string arg)
{
    return "exp(" + arg + ")";
}

double gauss(double a)
{
    return std::exp(-a * a);
}
double dgauss(double a)
{
    return -2. * a * std::exp(-a * a);
}
double ddgauss(double a)
{
    return (4 * a * a - 2.) * std::exp(-a * a);
}
std::string pgauss(std::string arg)
{
    return "exp(-" + arg + "*" + arg + ")";
}
} // namespace dsyre
