#pragma once
#include "matrix.hpp"
#include "operations.hpp"
#include <map>


namespace linear_algebra::decompositions {
    template<typename Field>
    std::map<char, Matrix<Field>> LU(const Matrix<Field>& A) {
        if (A.get_cols() != A.get_rows()) {
            throw std::runtime_error("Can't compute LU decomposition for non-square matrix");
        }
        size_t n = A.get_rows();
        Matrix<Field> L(n), U(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (j < i) {
                    U[i][j] = 0;
                } else {
                    U[i][j] = 0;
                    for (size_t k = 0; k < i; ++k) {
                        U[i][j] += L[i][k] * U[k][j];
                    }
                    U[i][j] = A[i][j] - U[i][j];
                }
            }
            for (size_t j = 0; j < n; ++j) {
                if (j < i) {
                    L[j][i] = 0;
                } else if (j == i) {
                    L[j][i] = 1;
                } else {
                    L[j][i] = 0;
                    for (size_t k = 0; k < i; ++k) {
                        L[j][i] += L[j][k] * U[k][i];
                    }
                    L[j][i] = (A[j][i] - L[j][i]) / U[i][i];
                }
            }
        }
        std::map<char, Matrix<Field>> result;
        result['L'] = L;
        result['U'] = U;
        return result;
    }

    template<typename Field>
    std::map<char, Matrix<Field>> QR(const Matrix<Field>& A) {
        Matrix<Field> Q(A);
        Q.transpose();
        size_t n = Q.get_rows();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                Q[i] -= dot_product(Q[i], Q[j]) * Q[j];
            }
            Q[i] /= norm2(Q[i]);
        }
        Matrix<Field> R = Q * A;
        std::map<char, Matrix<Field>> result;
        result['Q'] = Q.transpose();
        result['R'] = R;
        return result;
    }

    template<typename Field>
    std::map<char, Matrix<Field>> SVD(const Matrix<Field>& A) {
        size_t n = A.get_rows();
        size_t m = A.get_cols();
        Matrix<Field> U(n, n), D(n, m), V(m, m);
        auto tmp = eigen_values_and_vectors_for_normal_matrix(A.T() * A);
        for (size_t i = 0; i < m; ++i) {
            V[i] = tmp[i].second;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    D[j][i] = 0;
                } else {
                    D[j][i] = std::sqrt(tmp[i].first);
                }
            }
        }
        V.transpose();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; i < n; ++i) {
                U[i][j] = 0;
            }
        }
        for (size_t i = 0; (i < m) && (i < n); ++i) {
            if (D[i][i] != 0) {
                U[i] = (A * V[i]) / D[i][i];
            }
        }
        V.transpose();
        U.transpose();
        std::map<char, Matrix<Field>> result;
        result['U'] = U;
        result['D'] = D;
        result['V'] = V;
        return result;
    }
} // decompositions
