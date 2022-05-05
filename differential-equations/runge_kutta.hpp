#pragma once
#include <vector>
#include <utility>
#include <functional>

namespace differential_equations {
    template<typename F>
    std::vector<std::pair<double, double>> runge_kutta_for_eq_second_order(F func,
                                                                           double segment_len,
                                                                           uint32_t dots_count,
                                                                           double x0,
                                                                           double y0) {
        double h = segment_len / dots_count;
        std::vector<std::pair<double, double>> result = {std::make_pair(x0, y0)};
        for (uint32_t i = 0; i < dots_count; ++i) {
            y0 += (func(x0, y0) + func(x0 + h, y0 + h * func(x0, y0))) * h / 2;
            x0 += h;
            result.emplace_back(std::make_pair(x0, y0));
        }
        return result;
    }


    template<typename F>
    std::vector<std::pair<double, double>> runge_kutta_for_eq_fourth_order(F func,
                                                                           double segment_len,
                                                                           uint32_t dots_count,
                                                                           double x0,
                                                                           double y0) {
        double h = segment_len / dots_count;
        std::vector<std::pair<double, double>> result = {std::make_pair(x0, y0)};
        for (uint32_t i = 0; i < dots_count; ++i) {
            auto c1 = func(x0, y0);
            auto c2 = func(x0 + h / 2, y0 + h / 2 * c1);
            auto c3 = func(x0 + h / 2, y0 + h / 2 * c2);
            auto c4 = func(x0 + h, y0 + h * c3);
            y0 += (c1 + 2 * c2 + 2 * c3 + c4) * h / 6;
            x0 += h;
            result.emplace_back(std::make_pair(x0, y0));
        }
        return result;
    }
} // differential_equations
