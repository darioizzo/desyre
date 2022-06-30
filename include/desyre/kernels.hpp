// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the heyoka library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DESYRE_KERNELS_HPP
#define DESYRE_KERNELS_HPP

namespace desyre
{
using kernel_f_ptr = double (*)(double, double);

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


// Global array of function pointers
kernel_f_ptr kernel_list[] = {cos, sin, exp};
kernel_f_ptr dkernel_list[] = {dcos, dsin, dexp};
kernel_f_ptr ddkernel_list[] = {ddcos, ddsin, ddexp};

} // namespace desyre

#endif