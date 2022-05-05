#pragma once
#include <random>
#include <vector>
#include <utility>
#include <iostream>

namespace integration {
    // Abstract class for area in multidimensional space
    class Area {
     public:
        // checks is given dot in the area
        virtual bool operator()(const std::vector<double>&) const {
            return true;
        }
        // the limits of the area for each axis (can be larger than the real one)
        virtual std::vector<std::pair<double, double>> get_limits() const {
            return {};
        }

        // volume of the area, if it's difficult to calculate it, return 0
        // In this case it will be calculated with Monte-Carlo method
        virtual double volume() const {
            return 0;
        }

        virtual ~Area() = default;
    };

    // Simple one dimensional Monte-Carlo method with uniform distribution
    template<typename F>
    double monte_carlo_integral(F func, double from, double to, uint32_t n) {
        double result = 0;
        auto dis = std::uniform_real_distribution(from, to);
        std::random_device rd;
        std::mt19937 gen(rd());
        for (uint32_t i = 0; i < n; ++i) {
            result += func(dis(gen));
        }
        return (to - from) * result / n;
    }

    // Simple  multidimensional Monte-Carlo method with uniform distribution
    template<typename F>
    double monte_carlo_integral(F func, const Area& area, uint32_t d, uint32_t n) {
        auto limits = area.get_limits();
        double jacobian = 1;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::vector<std::uniform_real_distribution<double>> dis;
        for (auto& limit : limits) {
            jacobian *= limit.second - limit.first;
            dis.emplace_back(std::uniform_real_distribution(0.0, 1.0));
        }
        double result = 0;
        std::vector<double> dot;
        dot.resize(d);
        uint32_t count = 0;
        for (uint32_t i = 0; i < n; ++i) {
            for (int j = 0; j < d; ++j) {
                dot[j] = limits[j].first + (limits[j].second - limits[j].first) * dis[j](gen);
            }
            if (area(dot)) {
                result += func(dot);
                ++count;
            }
        }
        if (area.volume() != 0) {
            return result / count * area.volume();
        }
        return jacobian / n * result;
    }
}  // namespace integration
