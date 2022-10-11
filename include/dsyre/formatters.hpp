// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_FORMATTERS_HPP
#define DSYRE_FORMATTERS_HPP

#include <fmt/format.h>
#include <symengine/expression.h>

template <>
struct fmt::formatter<SymEngine::Expression> {
    // Presentation format: 'v' - vanilla, 'e' - expanded.
    char presentation = 'v';

    // Parses format specifications of the form ['s' | 'e'].
    constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin())
    {
        // [ctx.begin(), ctx.end()) is a character range that contains a part of
        // the format string starting from the format specifications to be parsed,
        // e.g. in
        //
        //   fmt::format("{:s} - simplified expression", expression);
        //
        // the range will contain "s} - simplified expression". The formatter should
        // parse specifiers until '}' or the end of the range. The formatter parses
        // the 's' specifier and returns an iterator pointing to '}'.

        // Please also note that this character range may be empty, in case of
        // the "{}" format string, so therefore you should check ctx.begin()
        // for equality with ctx.end().

        // Parse the presentation format and stores it in the formatter:
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'v' || *it == 'e')) {
            presentation = *it++;
        }
        // Check if reached the end of the range:
        if (it != end && *it != '}') {
            throw format_error("invalid format");
        }
        // Return an iterator past the end of the parsed range:
        return it;
    }

    // Formats the expression ex using the parsed format specification (presentation)
    // stored in this formatter.
    template <typename FormatContext>
    auto format(const SymEngine::Expression &ex, FormatContext &ctx) const -> decltype(ctx.out())
    {
        std::ostringstream oss;
        presentation == 'v' ? oss << ex : oss << SymEngine::expand(ex);

        //  ctx.out() is an output iterator to write to.
        return fmt::format_to(ctx.out(), "{}", oss.str());
    }
};

#endif