// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the heyoka library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <dsyre/kernels.hpp>

namespace dsyre
{
// Functions and derivatives
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
} // namespace dsyre
