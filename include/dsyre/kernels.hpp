// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_KERNELS_HPP
#define DSYRE_KERNELS_HPP

#include <string>

namespace dsyre
{
// Functions and derivatives

double inv(double a);
double dinv(double a);
double ddinv(double a);
std::string pinv(std::string);

double cos(double a);
double dcos(double a);
double ddcos(double a);
std::string pcos(std::string);

double sin(double a);
double dsin(double a);
double ddsin(double a);
std::string psin(std::string);

double exp(double a);
double dexp(double a);
double ddexp(double a);
std::string pexp(std::string);

} // namespace dsyre

#endif