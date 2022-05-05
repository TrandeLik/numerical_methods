#pragma once
#include <cmath>

namespace root_finding {
    template<typename F>
    double secant_method(F func, double from, double to, double eps) {
        while (std::abs(to - from) >= eps) {
            double tmp = to;
            to = from - func(from) * (to - from) / (func(to) - func(from));
            from = tmp;
        }
        return to;
    }
} // root_finding
