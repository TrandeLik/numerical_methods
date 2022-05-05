#pragma once
#include <vector>
#include <functional>
#include <utility>

namespace differential_equations {
    template<typename F>
    std::vector<std::pair<double, double>> run_throw_for_boundary_value_problem(F p, F q, F f,
                                                                                double from, double to,
                                                                                double alpha1, double betta1,
                                                                                double gamma1,
                                                                                double alpha2, double betta2,
                                                                                double gamma2,
                                                                                uint32_t n) {
        std::vector<double> x(n + 1, 0);
        std::vector<double> y(n + 1, 0);
        std::vector<double> alpha(n + 1, 0);
        std::vector<double> betta(n + 1, 0);
        x[0] = from;
        x[n] = to;
        double h = (to - from) / n;
        alpha[1] = -betta1 / (alpha1 * h - betta1);
        betta[1] = gamma1 / (alpha1 - betta1 / h);
        for (uint32_t i = 1; i < n; ++i) {
            x[i] = x[i - 1] + h;
            auto A = 1 / (h * h) + p(x[i]) / (2 * h);
            auto B = 1 / (h * h) - p(x[i]) / (2 * h);
            auto C = -2 / (h * h) + q(x[i]);
            alpha[i + 1] = -A / (B * alpha[i] + C);
            betta[i + 1] = (f(x[i]) - B * betta[i]) / (B * alpha[i] + C);
        }
        y[n] = (betta2 * betta[n] + gamma2 * h) / (betta2 * (1 - alpha[n]) + alpha2 * h);
        for (uint32_t i = n; i != 0; --i) {
            y[i - 1] = y[i] * alpha[i] + betta[i];
        }
        std::vector<std::pair<double, double>> result(n + 1);
        for (int i = 0; i < n + 1; ++i) {
            result[i] = std::make_pair(x[i], y[i]);
        }
        return result;
    }
} // differential_equations
