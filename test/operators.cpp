// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>
#include <stdexcept>

#include <dsyre/expression.hpp>
#include <dsyre/kernels.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <dsyre/symengine.hpp>

#include "catch.hpp"

using namespace dsyre;

TEST_CASE("simple_check")
{
    expression ex(1u, 1u, {"sum", "diff", "mul", "div", "inv", "cos", "sin", "exp", "gauss"});
    std::vector<unsigned> geno;
    std::vector<double> vars{1.23};
    std::vector<double> cons{-0.12};
    std::vector<double> retval, dretval, ddretval;

    // +
    geno = std::vector<unsigned>{0, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == vars[0] + cons[0]);
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == 1.);
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == 0.);

    // -
    geno = std::vector<unsigned>{1, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == vars[0] - cons[0]);
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == 1.);
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == 0.);

    // *
    geno = std::vector<unsigned>{2, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == vars[0] * cons[0]);
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == cons[0]);
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == 0.);

    // /
    geno = std::vector<unsigned>{3, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == vars[0] / cons[0]);
    ex.dphenotype(dretval, geno, retval, 1); // c
    REQUIRE(dretval[2] == -vars[0] / cons[0] / cons[0]);
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == 2. * vars[0] / cons[0] / cons[0] / cons[0]);

    // inv
    geno = std::vector<unsigned>{4, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == 1 / vars[0]);
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == -1. / vars[0] / vars[0]);
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == 2. * 1. / vars[0] / vars[0] / vars[0]);

    // cos
    geno = std::vector<unsigned>{5, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == std::cos(vars[0]));
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == -std::sin(vars[0]));
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == -std::cos(vars[0]));

    // sin
    geno = std::vector<unsigned>{6, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == std::sin(vars[0]));
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == std::cos(vars[0]));
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == -std::sin(vars[0]));

    // exp
    geno = std::vector<unsigned>{7, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == std::exp(vars[0]));
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == std::exp(vars[0]));
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == std::exp(vars[0]));

    // gauss
    geno = std::vector<unsigned>{8, 0, 1};
    ex.phenotype(retval, geno, vars, cons);
    REQUIRE(retval[2] == std::exp(-vars[0] * vars[0]));
    ex.dphenotype(dretval, geno, retval, 0);
    REQUIRE(dretval[2] == -2. * vars[0] * std::exp(-vars[0] * vars[0]));
    ex.ddphenotype(ddretval, geno, retval, dretval, dretval);
    REQUIRE(ddretval[2] == -2 * std::exp(-vars[0] * vars[0]) + 4 * vars[0] * vars[0] * std::exp(-vars[0] * vars[0]));
}