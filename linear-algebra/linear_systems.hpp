#pragma once
#include "matrix.hpp"
#include "operations.hpp"


namespace linear_algebra::linear_systems {
    template<typename Field>
    Vector<Field> gauss_method(Matrix<Field> A, Vector<Field> b) {
        if (A.get_cols() != A.get_rows()) {
            throw std::runtime_error("Can't solve system with non square matrix");
        }
        if (A.get_rows() != b.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t n = A.get_cols();
        for (size_t i = 0; i < n; ++i) {
            size_t max_index = i;
            for (size_t j = i; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[max_index][i])) {
                    max_index = j;
                }
            }
            if (A[max_index][i] == 0) {
                throw std::runtime_error("The matrix of the system is degenerate");
            }
            swap_rows(A, i, max_index);
            std::swap(b[i], b[max_index]);
            for (size_t j = i + 1; j < n; ++j) {
                auto ratio = A[j][i] / A[i][i];
                A[j] -= ratio * A[i];
                b[j] -= ratio * b[i];
            }
        }
        Vector<Field> x(n);
        x[n - 1] = b[n - 1] / A[n - 1][n - 1];
        for (size_t i = n - 2; i >= 0; --i) {
            x[i] = b[i];
            for (size_t j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] = x[i] / A[i][i];
            if (i == 0) {
                break;
            }
        }
        return x;
    }

    template<typename Field>
    Vector<Field>
    upper_relaxation(const Matrix<Field> &A, const Vector<Field> &b, double w = 1.7, double eps = 0.00000001,
                     uint32_t max_iter = 100000) {
        if (A.get_cols() != A.get_rows()) {
            throw std::runtime_error("Can't solve system with non square matrix");
        }
        if (A.get_rows() != b.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t n = A.get_rows();
        Vector<Field> x_prev(n);
        Vector<Field> x_cur(n);
        for (size_t i = 0; i < n; ++i) {
            x_cur[i] = 0;
        }
        uint32_t current_iter = 0;
        do {
            x_prev = x_cur;
            for (size_t i = 0; i < n; ++i) {
                Field c = 0;
                for (size_t j = 0; j < i; ++j) {
                    c += A[i][j] * x_cur[j];
                }
                for (size_t j = i; j < n; ++j) {
                    c += A[i][j] * x_prev[j];
                }
                x_cur[i] = x_prev[i] + w * (b[i] - c) / A[i][i];
            }
            ++current_iter;
            if (current_iter > max_iter) {
                throw std::runtime_error("Iteration limit exceeded");
            }
        } while (norm2(x_prev - x_cur) >= eps / 2);
        return x_cur;
    }

    template<typename Field>
    Vector<Field> steepest_descent(const Matrix<Field> &A, const Vector<Field> &b, double eps = 0.00000001,
                                   uint32_t max_iter = 100000) {
        if (A.get_rows() != b.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t m = A.get_rows();
        size_t n = A.get_cols();
        Vector<Field> x(n);
        for (size_t i = 0; i < n; ++i) {
            x[i] = 0;
        }
        Vector<Field> r = A * x - b;
        size_t current_iter = 0;
        while (norm2(r) >= eps) {
            Field alpha = dot_product(r, r) / dot_product(A * r, r);
            x -= alpha * r;
            r = A * x - b;
            ++current_iter;
            if (current_iter > max_iter) {
                throw std::runtime_error("Iteration limit exceeded");
            }
        }
        return x;
    }

    template<typename Field>
    Vector<Field> minimal_residual(const Matrix<Field> &A, const Vector<Field> &b, double eps = 0.00000001,
                                   uint32_t max_iter = 100000) {
        if (A.get_rows() != b.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t m = A.get_rows();
        size_t n = A.get_cols();
        Vector<Field> x(n);
        for (size_t i = 0; i < n; ++i) {
            x[i] = 0;
        }
        Vector<Field> r = A * x - b;
        size_t current_iter = 0;
        while (norm2(r) >= eps) {
            Field alpha = dot_product(r, r) / dot_product(A * r, A * r);
            x -= alpha * r;
            r = A * x - b;
            ++current_iter;
            if (current_iter > max_iter) {
                throw std::runtime_error("Iteration limit exceeded");
            }
        }
        return x;
    }

    template<typename Field>
    Vector<Field> conjugate_gradient(const Matrix<Field> &A, const Vector<Field> &b, double eps = 0.00000001) {
        if (A.get_rows() != b.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t n = A.get_cols();
        Vector<Field> x(n);
        for (size_t i = 0; i < n; ++i) {
            x[i] = 0;
        }
        Vector<Field> r = b - A * x;
        Vector<Field> p(r);
        while (norm2(r) >= eps) {
            Field alpha = dot_product(p, r) / dot_product(p, A * p);
            x += alpha * p;
            r -= alpha * A * p;
            Field betta = -dot_product(p, A * r) / dot_product(p, A * p);
            p = r + betta * p;
        }
        return x;
    }
} // linear_systems
