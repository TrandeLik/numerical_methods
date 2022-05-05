#pragma once
#include <cmath>

namespace root_finding {
    template<typename F, typename G>
    double newtons_method(F func, G derivative, double from, double to, double eps) {
        double x = (from + to) / 2;
        while (std::abs(func(x)) >= eps / 2) {
            x -= func(x) / derivative(x);
        }
        return x;
    }

    template<typename F>
    double _derivative(F func, double x, double h) {
        return (func(x + h) - func(x - h)) / (2 * h);
    }

    template<typename F>
    double newtons_method(F func, double from, double to, double eps) {
        double x = (from + to) / 2;
        while (std::abs(func(x)) >= eps / 2) {
            x -= func(x) / _derivative(func, x, eps);
        }
        return x;
    }
} // root_finding
