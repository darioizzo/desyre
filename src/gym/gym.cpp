// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <vector>

#include <dsyre/detail/data_generation.hpp>
#include <dsyre/detail/visibility.hpp>
#include <dsyre/gym/gym.hpp>
#include <dsyre/gym/misc_data.hpp>
#include <dsyre/gym/nist_data.hpp>

namespace dsyre
{
namespace gym
{
DSYRE_DLL_PUBLIC void generate_P0(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::koza_quintic, -3., 3., 10);
}
DSYRE_DLL_PUBLIC void generate_P1(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P1, 1., 3., 10);
}
DSYRE_DLL_PUBLIC void generate_P2(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P2, 0.1, 5., 10);
}
DSYRE_DLL_PUBLIC void generate_P3(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P3, -0.9, 1, 10);
}
DSYRE_DLL_PUBLIC void generate_P4(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P4, -1, 1, 10);
}
DSYRE_DLL_PUBLIC void generate_P5(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P5, 1., 3., 10);
}
DSYRE_DLL_PUBLIC void generate_P6(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P6, -2.1, 1., 10);
}
DSYRE_DLL_PUBLIC void generate_P7(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P7, -1, 1, 10);
}
DSYRE_DLL_PUBLIC void generate_P8(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_normal_data(points, labels, detail::P8, 100, 5);
}
DSYRE_DLL_PUBLIC void generate_P9(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.3, 4);
    for (unsigned i = 0; i < 100u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::kotanchek(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P10(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::salutowicz, 0.05, 10., 100);
}
DSYRE_DLL_PUBLIC void generate_P11(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 10);
    for (unsigned i = 0; i < 601u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::salutowicz2d(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P12(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 1024u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt), dist(mt), dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::uball5d(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P13(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 2);
    std::uniform_real_distribution<double> dist1(1, 2);
    for (unsigned i = 0; i < 300u; ++i) {
        std::vector<double> point = {dist(mt), dist1(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::ratpol3d(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P14(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.1, 5.9);
    for (unsigned i = 0; i < 30u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::sinecosine(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P15(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 300u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::ripple(point)});
    }
}
DSYRE_DLL_PUBLIC void generate_P16(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 50u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::ratpol2d(point)});
    }
}

DSYRE_DLL_PUBLIC void generate_P17(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = chwirut1_points;
    labels = chwirut1_labels;
}

DSYRE_DLL_PUBLIC void generate_P18(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = chwirut2_points;
    labels = chwirut2_labels;
}

DSYRE_DLL_PUBLIC void generate_P19(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = daniel_wood_points;
    labels = daniel_wood_labels;
}

DSYRE_DLL_PUBLIC void generate_P20(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = gauss1_points;
    labels = gauss1_labels;
}

DSYRE_DLL_PUBLIC void generate_P21(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = kirby2_points;
    labels = kirby2_labels;
}

DSYRE_DLL_PUBLIC void generate_P22(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = lanczos2_points;
    labels = lanczos2_labels;
}

DSYRE_DLL_PUBLIC void generate_P23(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = misra1b_points;
    labels = misra1b_labels;
}

DSYRE_DLL_PUBLIC void generate_P24(std::vector<std::vector<double>> &points, std::vector<double> &labels)
{
    points.clear();
    labels.clear();
    points = luca1_points;
    labels = luca1_labels;
}

} // namespace gym
} // namespace dsyre