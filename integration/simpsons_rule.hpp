#pragma once
#include <functional>

namespace integration {
    template<typename F>
    double simpsons_integral(F func, double from, double to, double eps) {
        double previous, current = 0;
        uint32_t n = 1;
        do {
            n *= 2;
            previous = current;
            current = 0;
            double h = (to - from) / n;
            double left = from, right;
            for (uint32_t i = 0; i < n; ++i) {
                right = left + h;
                current += (right - left) / 6 * (func(left) + 4 * func((left + right) / 2) + func(right));
                left = right;
            }
        } while (std::abs(previous - current) >= eps / 2);
        return current;
    }
}  // namespace integration
