#pragma once
#include <cmath>
#include <utility>
#include <random>
#include <vector>
#include "matrix.hpp"

namespace linear_algebra {
    template <typename Field>
    Field norm2(const Vector<Field>& v) {
        Field result = 0;
        size_t n = v.get_size();
        for (int i = 0; i < n; ++i) {
            result += v[i] * v[i];
        }
        return std::sqrt(result);
    }

    template <typename Field>
    Field dot_product(const Vector<Field>& v1, const Vector<Field>& v2) {
        if (v1.get_size() != v2.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t n = v1.get_size();
        Field result = 0;
        for (size_t i = 0; i < n; ++i) {
            result += v1[i] * v2[i];
        }
        return result;
    }

    template <typename Field>
    std::pair<Field, Vector<Field>> greatest_eigen_value_and_vector(const Matrix<Field>& A, double eps=0.00001, uint32_t max_iter=100000) {
        if (A.get_cols() != A.get_rows()) {
            throw std::runtime_error("Can't calculate eigen value for non-square matrix");
        }
        auto dis = std::uniform_real_distribution(0.0, 1.1);
        std::random_device rd;
        std::mt19937 gen(rd());
        size_t n = A.get_rows();
        uint32_t current_iter = 0;
        Vector<Field> x(n);
        bool flag = true;
        do {
            for (size_t i = 0; i < n; ++i) {
                x[i] = dis(gen);
            }
            Vector<Field> y = A * x;
            for (size_t i = 0; i < n - 1; ++i) {
                if (!(x[i] == 0 && y[i] == 0) &&
                    !(x[i + 1] == 0 && y[i + 1] == 0) &&
                    (std::abs(y[i] / x[i] - y[i + 1] / y[i]) > 0.0001)) {
                    flag = false;
                }
            }
        } while (flag);
        Field eigen_prev = 0, eigen_cur = 0;
        do {
            eigen_prev = eigen_cur;
            eigen_cur = dot_product(A * x, x) / dot_product(x, x);
            x = A * x;
            x /= norm2(x);
            ++current_iter;
            if (current_iter > max_iter) {
                throw std::runtime_error("Iteration limit exceeded");
            }
        } while (std::abs(eigen_cur - eigen_prev) >= eps / 2);
        return std::make_pair(eigen_cur, x);
    }

    template <typename Field>
    Field greatest_eigen_value(const Matrix<Field>& A, double eps=0.00001, uint32_t max_iter=100000) {
        return greatest_eigen_value_and_vector(A, eps, max_iter).first;
    }

    template <typename Field>
    Field singular_norm(const Matrix<Field>& A) {
        return std::sqrt(greatest_eigen_value(A * A.T()));
    }

    template <typename Field>
    Field conditional_number(const Matrix<Field>& A) {
        return singular_norm(A) * singular_norm(A.inv());
    }

    template <typename Field>
    std::vector<std::pair<Field, Vector<Field>>> eigen_values_and_vectors_for_normal_matrix(Matrix<Field> A, double eps=0.00001, uint32_t max_iter=100000) {
        if (A.get_rows() != A.get_cols()) {
            throw std::runtime_error("Can't compute eigen values and vectors for non-square matrix");
        }
        std::vector<std::pair<Field, Vector<Field>>> result;
        size_t n = A.get_rows();
        for (size_t i = 0; i < n; ++i) {
            result.template emplace_back(greatest_eigen_value_and_vector(A, eps, max_iter));
            auto x = static_cast<Matrix<Field>>(result[i].second);
            A -= result[i].first * x * x.T();
        }
        return result;
    }
}
