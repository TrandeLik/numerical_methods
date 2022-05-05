#pragma once
#include <iostream>
#include <cstring>
#include <cstdint>
#include <utility>

namespace linear_algebra {
    template<typename Field = double>
    class Matrix;

    template<typename Field = double>
    class Vector {
    private:
        size_t size;
        Field *data_;
    public:
        Vector() : size(0), data_(nullptr) {}

        explicit Vector(size_t n);

        Vector(const Vector<Field> &r);

        size_t get_size() const;

        Field &operator[](size_t i);

        Field operator[](size_t i) const;

        Vector<Field> &operator=(Vector<Field> r);

        Vector<Field> &operator*=(Field x);

        Vector<Field> &operator/=(Field x);

        Vector<Field> &operator+=(const Vector<Field> &r);

        Vector<Field> &operator-=(const Vector<Field> &r);

        bool operator==(const Vector<Field> &r2) const;

        bool operator!=(const Vector<Field> &r2) const;

        explicit operator Matrix<Field>() const;

        template<typename F>
        friend std::ostream &operator<<(std::ostream &out, const Vector<F> &r);

        ~Vector();
    };

    template<typename Field>
    class Matrix {
        size_t rows_count;
        size_t cols_count;
        Vector<Field> *rows_;

    public:
        Matrix() : rows_count(0), cols_count(0), rows_(nullptr) {}

        Matrix(size_t n, size_t m);

        explicit Matrix(size_t n);

        Matrix(const Matrix &m);

        size_t get_rows() const;

        size_t get_cols() const;

        Vector<Field> &operator[](size_t i);

        Vector<Field> operator[](size_t i) const;

        Matrix<Field> &operator=(Matrix<Field> m);

        Matrix<Field> &operator*=(Field x);

        Matrix<Field> &operator/=(Field x);

        Matrix<Field> &operator+=(const Matrix<Field> &m);

        Matrix<Field> &operator-=(const Matrix<Field> &m);

        Matrix<Field> &operator*=(const Matrix<Field> &m);

        template<typename F>
        friend Matrix<F> operator*(const Matrix<F> &m1, const Matrix<F> &m2);

        bool operator==(const Matrix<Field> &m) const;

        bool operator!=(const Matrix<Field> &m) const;

        explicit operator Vector<Field>() const;

        Matrix<Field> T() const;

        Matrix<Field> &transpose();

        Field det() const;

        Matrix<Field> inv() const;

        Field trace() const;

        size_t rank() const;

        template<typename F>
        friend std::ostream &operator<<(std::ostream &out, const Matrix<F> &m);

        ~Matrix();
    };

    template<typename Field>
    Vector<Field>::Vector(size_t n) {
        size = n;
        data_ = new Field[size];
    }

    template<typename Field>
    Vector<Field>::Vector(const Vector<Field> &r) {
        size = r.size;
        data_ = new Field[size];
        if (r.data_ != nullptr) {
            memcpy(data_, r.data_, size * sizeof(*data_));
        }
    }

    template<typename Field>
    Vector<Field>::~Vector() {
        delete[] data_;
    }

    template<typename Field>
    Vector<Field> &Vector<Field>::operator=(Vector<Field> r) {
        std::swap(size, r.size);
        std::swap(data_, r.data_);
        return *this;
    }

    template<typename Field>
    Field &Vector<Field>::operator[](size_t i) {
        if (i >= size) {
            throw std::runtime_error("Bad index");
        }
        return data_[i];
    }

    template<typename Field>
    Field Vector<Field>::operator[](size_t i) const {
        if (i >= size) {
            throw std::runtime_error("Bad index");
        }
        return data_[i];
    }

    template<typename Field>
    size_t Vector<Field>::get_size() const {
        return size;
    }

    template<typename Field>
    Vector<Field> &Vector<Field>::operator*=(Field x) {
        for (size_t i = 0; i < size; ++i) {
            data_[i] *= x;
        }
        return *this;
    }

    template<typename Field>
    Vector<Field> &Vector<Field>::operator/=(Field x) {
        for (size_t i = 0; i < size; ++i) {
            data_[i] /= x;
        }
        return *this;
    }

    template<typename Field>
    Vector<Field> operator*(const Vector<Field> &r1, Field x) {
        Vector tmp(r1);
        tmp *= x;
        return tmp;
    }

    template<typename Field>
    Vector<Field> operator*(Field x, const Vector<Field> &r1) {
        Vector tmp(r1);
        tmp *= x;
        return tmp;
    }

    template<typename Field>
    Vector<Field> operator/(const Vector<Field> &r1, Field x) {
        Vector tmp(r1);
        tmp /= x;
        return tmp;
    }

    template<typename Field>
    Vector<Field> &Vector<Field>::operator+=(const Vector &r) {
        if (size != r.size) {
            throw std::runtime_error("Mismatched sizes");
        }
        for (size_t i = 0; i < size; ++i) {
            data_[i] += r.data_[i];
        }
        return *this;
    }

    template<typename Field>
    Vector<Field> &Vector<Field>::operator-=(const Vector &r) {
        if (size != r.size) {
            throw std::runtime_error("Mismatched sizes");
        }
        for (size_t i = 0; i < size; ++i) {
            data_[i] -= r.data_[i];
        }
        return *this;
    }

    template<typename Field>
    Vector<Field> operator+(const Vector<Field> &r1, const Vector<Field> &r2) {
        Vector tmp(r1);
        tmp += r2;
        return tmp;
    }

    template<typename Field>
    Vector<Field> operator-(const Vector<Field> &r1, const Vector<Field> &r2) {
        Vector tmp(r1);
        tmp -= r2;
        return tmp;
    }

    template<typename Field>
    bool Vector<Field>::operator==(const Vector<Field> &r2) const {
        if (size != r2.size) {
            return false;
        }
        for (size_t i = 0; i < size; ++i) {
            if (data_[i] != r2.data_[i]) {
                return false;
            }
        }
        return true;
    }

    template<typename Field>
    bool Vector<Field>::operator!=(const Vector<Field> &r2) const {
        return !(*this == r2);
    }

    template<typename Field>
    Vector<Field>::operator Matrix<Field>() const {
        Matrix<Field> m(size, 1);
        for (size_t i = 0; i < size; ++i) {
            m[i][0] = data_[i];
        }
        return m;
    }

    template<typename Field>
    std::ostream &operator<<(std::ostream &out, const Vector<Field> &r) {
        for (size_t i = 0; i < r.size; ++i) {
            out << r.data_[i] << ' ';
        }
        return out;
    }

//===============MATRIX===============
    template<typename Field>
    Matrix<Field>::Matrix(size_t n, size_t m) {
        rows_count = n;
        cols_count = m;
        rows_ = reinterpret_cast<Vector<Field> *>(new char[rows_count * sizeof(Vector<Field>)]);
        for (size_t i = 0; i < rows_count; ++i) {
            new(rows_ + i) Vector(cols_count);
        }
    }

    template<typename Field>
    Matrix<Field>::Matrix(size_t n) {
        rows_count = n;
        cols_count = n;
        rows_ = reinterpret_cast<Vector<Field> *>(new char[rows_count * sizeof(Vector<Field>)]);
        for (size_t i = 0; i < rows_count; ++i) {
            new(rows_ + i) Vector(cols_count);
        }
    }

    template<typename Field>
    Matrix<Field>::Matrix(const Matrix<Field> &m) {
        rows_count = m.rows_count;
        cols_count = m.cols_count;
        rows_ = reinterpret_cast<Vector<Field> *>(new char[rows_count * sizeof(Vector<Field>)]);
        for (size_t i = 0; i < rows_count; ++i) {
            new(rows_ + i) Vector(m.rows_[i]);
        }
    }

    template<typename Field>
    Matrix<Field>::~Matrix() {
        if (rows_ != nullptr) {
            for (size_t i = 0; i < rows_count; ++i) {
                (rows_ + i)->~Vector();
            }
            delete[] reinterpret_cast<char *>(rows_);
        }
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator=(Matrix<Field> m) {
        std::swap(rows_count, m.rows_count);
        std::swap(cols_count, m.cols_count);
        std::swap(rows_, m.rows_);
        return *this;
    }

    template<typename Field>
    size_t Matrix<Field>::get_rows() const {
        return rows_count;
    }

    template<typename Field>
    size_t Matrix<Field>::get_cols() const {
        return cols_count;
    }

    template<typename Field>
    Vector<Field> &Matrix<Field>::operator[](size_t i) {
        if (i >= rows_count) {
            throw std::runtime_error("Bad index");
        }
        return rows_[i];
    }

    template<typename Field>
    Vector<Field> Matrix<Field>::operator[](size_t i) const {
        if (i >= rows_count) {
            throw std::runtime_error("Bad index");
        }
        return rows_[i];
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator*=(Field x) {
        for (size_t i = 0; i < rows_count; ++i) {
            rows_[i] *= x;
        }
        return *this;
    }

    template<typename Field>
    Matrix<Field> operator*(const Matrix<Field> &m1, Field x) {
        Matrix tmp(m1);
        tmp *= x;
        return tmp;
    }

    template<typename Field>
    Matrix<Field> operator*(Field x, const Matrix<Field> &m1) {
        Matrix tmp(m1);
        tmp *= x;
        return tmp;
    }

    template<typename Field>
    Matrix<Field> operator/(const Matrix<Field> &m1, Field x) {
        Matrix tmp(m1);
        tmp /= x;
        return tmp;
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator/=(Field x) {
        for (size_t i = 0; i < rows_count; ++i) {
            rows_[i] /= x;
        }
        return *this;
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator+=(const Matrix<Field> &m) {
        if (m.rows_count != rows_count || m.cols_count != cols_count) {
            throw std::runtime_error("Mismatched sizes");
        }
        for (size_t i = 0; i < rows_count; ++i) {
            rows_[i] += m.rows_[i];
        }
        return *this;
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator-=(const Matrix<Field> &m) {
        if (m.rows_count != rows_count || m.cols_count != cols_count) {
            throw std::runtime_error("Mismatched sizes");
        }
        for (size_t i = 0; i < rows_count; ++i) {
            rows_[i] -= m.rows_[i];
        }
        return *this;
    }

    template<typename Field>
    Matrix<Field> operator+(const Matrix<Field> &m1, const Matrix<Field> &m2) {
        Matrix tmp(m1);
        tmp += m2;
        return tmp;
    }

    template<typename Field>
    Matrix<Field> operator-(const Matrix<Field> &m1, const Matrix<Field> &m2) {
        Matrix tmp(m1);
        tmp -= m2;
        return tmp;
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::operator*=(const Matrix<Field> &m) {
        if (cols_count != m.rows_count) {
            throw std::runtime_error("Mismatched sizes");
        }
        Matrix<Field> copy(*this);
        (*this) = Matrix<Field>(rows_count, m.cols_count);
        for (size_t i = 0; i < rows_count; ++i) {
            for (size_t j = 0; j < m.cols_count; ++j) {
                rows_[i][j] = 0;
                for (size_t k = 0; k < copy.cols_count; ++k) {
                    rows_[i][j] += copy[i][k] * m[k][j];
                }
            }
        }
        return *this;
    }

    template<typename Field>
    Matrix<Field> operator*(const Matrix<Field> &m1, const Matrix<Field> &m2) {
        if (m1.cols_count != m2.rows_count) {
            throw std::runtime_error("Mismatched sizes");
        }
        Matrix<Field> result(m1.rows_count, m2.cols_count);
        for (size_t i = 0; i < m1.rows_count; ++i) {
            for (size_t j = 0; j < m2.cols_count; ++j) {
                result[i][j] = 0;
                for (size_t k = 0; k < m1.cols_count; ++k) {
                    result[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return result;
    }

    template<typename Field>
    bool Matrix<Field>::operator==(const Matrix<Field> &m) const {
        if (cols_count != m.cols_count || rows_count != m.rows_count) {
            return false;
        }
        for (size_t i = 0; i < m.rows_count; ++i) {
            if (rows_[i] != m.rows_[i]) {
                return false;
            }
        }
        return true;
    }

    template<typename Field>
    bool Matrix<Field>::operator!=(const Matrix<Field> &m) const {
        return !(*this == m);
    }

    template<typename Field>
    std::ostream &operator<<(std::ostream &out, const Matrix<Field> &m) {
        for (size_t i = 0; i < m.rows_count; ++i) {
            out << m.rows_[i] << '\n';
        }
        return out;
    }

    template<typename Field>
    Matrix<Field>::operator Vector<Field>() const {
        if (cols_count != 1) {
            throw std::runtime_error("Can't cast to vector");
        }
        Vector<Field> result(rows_count);
        for (size_t i = 0; i < rows_count; ++i) {
            result[i] = rows_[i][0];
        }
        return result;
    }

    template<typename Field>
    Matrix<Field> unite_matrix(size_t n) {
        Matrix result(n);
        for (size_t i = 0; i < n; ++i) {
            result[i][i] = 1;
        }
        return result;
    }

    template<typename Field>
    void swap_rows(Matrix<Field> &m, size_t i, size_t j) {
        std::swap(m[i], m[j]);
    }

    template<typename Field>
    void swap_cols(Matrix<Field> &m, size_t i, size_t j) {
        size_t limit = m.get_rows();
        for (size_t k = 0; k < limit; ++k) {
            std::swap(m[k][i], m[k][j]);
        }
    }

    template<typename Field>
    Matrix<Field> Matrix<Field>::T() const {
        Matrix<Field> m(cols_count, rows_count);
        for (size_t i = 0; i < cols_count; ++i) {
            for (size_t j = 0; j < rows_count; ++j) {
                m[i][j] = rows_[j][i];
            }
        }
        return m;
    }

    template<typename Field>
    Matrix<Field> &Matrix<Field>::transpose() {
        Matrix<Field> m(cols_count, rows_count);
        for (size_t i = 0; i < cols_count; ++i) {
            for (size_t j = 0; j < rows_count; ++j) {
                m[i][j] = rows_[j][i];
            }
        }
        (*this) = m;
        return (*this);
    }

    template<typename Field>
    Field Matrix<Field>::det() const {
        if (rows_count != cols_count) {
            throw std::runtime_error("Can't calculate determinant of non-square matrix");
        }
        Matrix<Field> copy(*this);
        Field result = 1;
        for (size_t i = 0; i < rows_count; ++i) {
            size_t max_row = i;
            for (size_t j = i; j < rows_count; ++j) {
                if (std::abs(copy[j][i]) > std::abs(copy[max_row][i])) {
                    max_row = i;
                }
            }
            if (copy[max_row][i] == 0) {
                return 0;
            }
            if (i != max_row) {
                swap_rows(copy, i, max_row);
                result *= -1;
            }
            result *= copy[i][i];
            for (size_t j = i + 1; j < rows_count; ++j) {
                auto ratio = copy[j][i] / copy[i][i];
                for (size_t k = i; k < rows_count; ++k) {
                    copy[j][k] -= ratio * copy[i][k];
                }
            }
        }
        return result;
    }

    template<typename Field>
    Field Matrix<Field>::trace() const {
        size_t limit = std::min(rows_count, cols_count);
        Field result = 0;
        for (size_t i = 0; i < limit; ++i) {
            result += rows_[i][i];
        }
        return result;
    }

    template<typename Field>
    size_t Matrix<Field>::rank() const {
        if (rows_count < cols_count) {
            return T().rank();
        }
        Matrix<Field> copy(*this);
        size_t result = cols_count;
        for (size_t i = 0; i < cols_count; ++i) {
            size_t max_row = i;
            for (size_t j = i; j < rows_count; ++j) {
                if (std::abs(copy[j][i]) > std::abs(copy[max_row][i])) {
                    max_row = i;
                }
            }
            if (copy[max_row][i] == 0) {
                --result;
            } else {
                swap_rows(copy, i, max_row);
                for (size_t j = i + 1; j < rows_count; ++j) {
                    auto ratio = copy[j][i] / copy[i][i];
                    for (size_t k = i; k < cols_count; ++k) {
                        copy[j][k] -= ratio * copy[i][k];
                    }
                }
            }
        }
        return result;
    }

    template<typename Field>
    Matrix<Field> Matrix<Field>::inv() const {
        if (rows_count != cols_count) {
            throw std::runtime_error("Can't inverse non-square matrix");
        }
        Matrix<Field> copy(*this);
        Matrix<Field> inverted = unite_matrix<Field>(rows_count);
        for (size_t i = 0; i < rows_count; ++i) {
            size_t max_row = i;
            for (size_t j = i; j < rows_count; ++j) {
                if (std::abs(copy[j][i]) > std::abs(copy[max_row][i])) {
                    max_row = i;
                }
            }
            if (copy[max_row][i] == 0) {
                throw std::runtime_error("The matrix is degenerate");
            }
            swap_rows(copy, i, max_row);
            swap_rows(inverted, i, max_row);
            for (size_t j = 0; j < rows_count; ++j) {
                if (j != i) {
                    auto ratio = copy[j][i] / copy[i][i];
                    copy[j] -= ratio * copy[i];
                    inverted[j] -= ratio * inverted[i];
                }
            }
            auto ratio = copy[i][i];
            copy[i] /= ratio;
            inverted[i] /= ratio;
        }
        return inverted;
    }

//===============MATVEC===============
    template<typename Field>
    Vector<Field> operator*(const Matrix<Field> &m, const Vector<Field> &v) {
        if (m.get_cols() != v.get_size()) {
            throw std::runtime_error("Mismatched sizes");
        }
        size_t size = m.get_rows();
        size_t inner_size = v.get_size();
        Vector<Field> result(size);
        for (size_t i = 0; i < size; ++i) {
            result[i] = 0;
            for (size_t j = 0; j < inner_size; ++j) {
                result[i] += m[i][j] * v[j];
            }
        }
        return result;
    }
} // linear-algebra
