// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include <dsyre/expression.hpp>

#include "catch.hpp"

using namespace dsyre;

TEST_CASE("construction")
{
    expression ex(1, 1, {0, 1, 2, 3, 4, 5});
    REQUIRE_NOTHROW(ex.random_constants(2.,3.));
}
