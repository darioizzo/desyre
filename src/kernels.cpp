// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the heyoka library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include<cmath>

#include<desyre/kernels.hpp>

namespace desyre
{
// Functions and derivatives
double sin(double a, double b)
{
    return std::sin(a);
}
double dsin(double a, double b)
{
    return std::cos(a);
}
double ddsin(double a, double b)
{
    return -std::sin(a);
}
double cos(double a, double b)
{
    return std::cos(a);
}
double dcos(double a, double b)
{
    return -std::sin(a);
}
double ddcos(double a, double b)
{
    return -std::cos(a);
}
double exp(double a, double b)
{
    return std::exp(a);
}
double dexp(double a, double b)
{
    return std::exp(a);
}
double ddexp(double a, double b)
{
    return std::exp(a);
}
} // namespace desyre
