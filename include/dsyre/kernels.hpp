// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the heyoka library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_KERNELS_HPP
#define DSYRE_KERNELS_HPP

namespace dsyre
{
// Functions and derivatives
double sin(double a, double b);
double dsin(double a, double b);
double ddsin(double a, double b);

double cos(double a, double b);
double dcos(double a, double b);
double ddcos(double a, double b);

double exp(double a, double b);
double dexp(double a, double b);
double ddexp(double a, double b);

} // namespace dsyre

#endif