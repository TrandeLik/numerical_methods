#include <gtest/gtest.h>
#include <sstream>
#include "matrix.hpp"

using namespace linear_algebra;

class TestMatrix: public testing::Test {
 protected:
    void SetUp() {
    }

    void TearDown() {
    }
};

template<typename Field = int32_t>
Matrix<Field> fill_rectangle_matrix(int n, int m) {
    Matrix<Field> res(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i][j] = (i + 1) * (j + 1);
        }
    }
    return res;
}

template<typename Field = int32_t>
Matrix<Field> fill_square_matrix(int n) {
    Matrix<Field> res(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i][j] = (i + 1) * (j + 1);
        }
    }
    return res;
}


TEST_F(TestMatrix, get_params) {
    Matrix<int32_t> m;
    ASSERT_EQ(0, m.get_rows());
    ASSERT_EQ(0, m.get_cols());
    m = fill_rectangle_matrix(3, 4);
    ASSERT_EQ(3, m.get_rows());
    ASSERT_EQ(4, m.get_cols());
    m = fill_square_matrix(3);
    ASSERT_EQ(3, m.get_rows());
    ASSERT_EQ(3, m.get_cols());
    m = fill_rectangle_matrix(1000, 1200);
    const Matrix a(m);
    ASSERT_EQ(1000, a.get_rows());
    ASSERT_EQ(1200, a.get_cols());
}

TEST_F(TestMatrix, correct_assignment) {
    Matrix m = fill_rectangle_matrix(2, 3);
    ASSERT_EQ(1, m[0][0]);
    ASSERT_EQ(2, m[0][1]);
    ASSERT_EQ(3, m[0][2]);
    ASSERT_EQ(2, m[1][0]);
    ASSERT_EQ(4, m[1][1]);
    ASSERT_EQ(6, m[1][2]);
    m = fill_rectangle_matrix(1000, 1200);
    const Matrix a(m);
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ((i + 1) * (j + 1), a[i][j]);
        }
    }
}

TEST_F(TestMatrix, incorrect_assignment) {
    Matrix m = fill_rectangle_matrix(3, 4);
    ASSERT_THROW(m[10][0], std::runtime_error);
    ASSERT_THROW(m[0][10], std::runtime_error);
    ASSERT_THROW(m[10][10], std::runtime_error);
    ASSERT_THROW(m[-1][3], std::runtime_error);
    ASSERT_THROW(m[2][-1], std::runtime_error);
    ASSERT_THROW(m[3][3], std::runtime_error);
    ASSERT_THROW(m[2][4], std::runtime_error);
    ASSERT_THROW(m[3][4], std::runtime_error);
    ASSERT_THROW(m[-1][-1], std::runtime_error);
    ASSERT_THROW(m[-1][4], std::runtime_error);
    ASSERT_THROW(m[3][-1], std::runtime_error);
}

TEST_F(TestMatrix, multiply_and_divide) {
    Matrix m = fill_rectangle_matrix(1000, 1200);
    m *= 2;
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ(2 * (i + 1) * (j + 1), m[i][j]);
        }
    }
    m /= 2;
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ((i + 1) * (j + 1), m[i][j]);
        }
    }
    m = m * 2;
    m *= 3;
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ(6 * (i + 1) * (j + 1), m[i][j]);
        }
    }
    m = m / 3;
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ(2 * (i + 1) * (j + 1), m[i][j]);
        }
    }
    m *= 0;
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1200; ++j) {
            ASSERT_EQ(0, m[i][j]);
        }
    }
}

TEST_F(TestMatrix, printing) {
    std::stringstream s;
    Matrix m = fill_square_matrix(1);
    s << m;
    ASSERT_EQ("1 \n", s.str());
    s.str("");
    m = fill_rectangle_matrix(2, 3);
    s << m;
    ASSERT_EQ("1 2 3 \n2 4 6 \n", s.str());
}

TEST_F(TestMatrix, comparison) {
    Matrix<int32_t> m;
    Matrix<int32_t> m2;
    ASSERT_TRUE(m == m2);
    ASSERT_FALSE(m != m2);
    m = fill_rectangle_matrix(3, 3);
    ASSERT_FALSE(m == m2);
    ASSERT_TRUE(m != m2);
    m2 = fill_square_matrix(1);
    ASSERT_FALSE(m == m2);
    ASSERT_TRUE(m != m2);
    m2 = fill_square_matrix(3);
    ASSERT_TRUE(m == m2);
    ASSERT_FALSE(m != m2);
    m[0][0] = 0;
    ASSERT_FALSE(m == m2);
    ASSERT_TRUE(m != m2);
}

TEST_F(TestMatrix, add) {
    Matrix m = fill_rectangle_matrix(3, 10);
    Matrix m2 = fill_rectangle_matrix(3, 10);
    Matrix m3 = fill_rectangle_matrix(3, 10);
    m3 *= 2;
    ASSERT_EQ(m3, m + m2);
    ASSERT_EQ(m3 - m, m2);
    m3 -= m;
    ASSERT_EQ(m3, m2);
    m = fill_rectangle_matrix(3, 3);
    ASSERT_THROW(m + m2, std::runtime_error);
    m = fill_rectangle_matrix(2, 10);
    ASSERT_THROW(m + m2, std::runtime_error);
    m = fill_rectangle_matrix(10, 10);
    ASSERT_THROW(m + m2, std::runtime_error);
    ASSERT_EQ(m2 + m3 * 0, m2);
}

TEST_F(TestMatrix, multiple_operations) {
    Matrix m1 = fill_rectangle_matrix(3, 3);
    const Matrix m2 = fill_square_matrix(3);
    Matrix m3 = m1 + m2;
    m3 *= 2;
    m1 *= 3;
    (m3 += m1) *= -1;
    m3[0][0] = 10000;
    std::stringstream s;
    s << m3;
    ASSERT_EQ("10000 -14 -21 \n-14 -28 -42 \n-21 -42 -63 \n", s.str());
}

TEST_F(TestMatrix, number) {
    Matrix m1 = fill_square_matrix(1);
    Matrix m2 = fill_square_matrix(1);
    m2[0][0] = 3;
    ASSERT_TRUE(m1 != m2);
    ASSERT_EQ(m2, m1 + m1 + m1);
    m1 *= 3;
    ASSERT_EQ(m2, m1);
    ASSERT_TRUE(m1 == m2);
    ASSERT_EQ(1, m1.get_cols());
    ASSERT_EQ(1, m1.get_rows());
}

TEST_F(TestMatrix, vector) {
    Matrix m1 = fill_rectangle_matrix(10, 1);
    Matrix m2 = fill_rectangle_matrix(10, 1);
    m2 *= 3;
    ASSERT_TRUE(m1 != m2);
    ASSERT_EQ(m2, m1 + m1 + m1);
    m1 *= 3;
    ASSERT_EQ(m2, m1);
    ASSERT_TRUE(m1 == m2);
    ASSERT_EQ(1, m1.get_cols());
    ASSERT_EQ(10, m1.get_rows());
}

TEST_F(TestMatrix, transposed_vector) {
    Matrix m1 = fill_rectangle_matrix(1, 10);
    Matrix m2 = fill_rectangle_matrix(1, 10);
    m2 *= 3;
    ASSERT_TRUE(m1 != m2);
    ASSERT_EQ(m2, m1 + m1 + m1);
    m1 *= 3;
    ASSERT_EQ(m2, m1);
    ASSERT_TRUE(m1 == m2);
    ASSERT_EQ(10, m1.get_cols());
    ASSERT_EQ(1, m1.get_rows());
}

TEST_F(TestMatrix, matrix_transpose) {
    Matrix m = fill_rectangle_matrix(5, 3);
    Matrix m2 = m.T();
    for (size_t i = 0; i < m2.get_rows(); ++i) {
        for (size_t j = 0; j < m2.get_cols(); ++j) {
            ASSERT_EQ((j + 1) * (i + 1), m2[i][j]);
        }
    }
    m.transpose();
    ASSERT_TRUE(m == m2);
}

TEST_F(TestMatrix, trace) {
    Matrix m = fill_rectangle_matrix(5, 3);
    ASSERT_TRUE(std::abs(14.0 - m.trace()) < 0.000001);
    m.transpose();
    ASSERT_TRUE(std::abs(14.0 - m.trace()) < 0.000001);
}

TEST_F(TestMatrix, determinant) {
    Matrix m(2, 2);
    m[0][0] = 2;
    m[0][1] = 1;
    m[1][0] = 1;
    m[1][1] = 4;
    ASSERT_TRUE(std::abs(7.0 - m.det()) < 0.000001);
    Matrix m1(3, 3);
    m1[0][0] = 2;
    m1[0][1] = 1;
    m1[1][0] = 1;
    m1[1][1] = 4;
    m1[2][2] = 5;
    ASSERT_TRUE(std::abs(35.0 - m1.det()) < 0.000001);
    Matrix m2(3, 2);
    ASSERT_THROW(m2.det(), std::runtime_error);
    m2 = fill_square_matrix<double>(10);
    ASSERT_TRUE(std::abs(m2.det()) < 0.000001);
}

TEST_F(TestMatrix, rank) {
    Matrix m(2, 2);
    m[0][0] = 2;
    m[0][1] = 1;
    m[1][0] = 1;
    m[1][1] = 4;
    ASSERT_EQ(2, m.rank());
    Matrix m1(3, 3);
    m1[0][0] = 2;
    m1[0][1] = 1;
    m1[1][0] = 1;
    m1[1][1] = 4;
    m1[2][2] = 5;
    ASSERT_EQ(3, m1.rank());
    Matrix m2 = fill_rectangle_matrix<double>(2, 3);
    ASSERT_EQ(1, m2.rank());
    m2 = fill_rectangle_matrix<double>(3, 2);
    ASSERT_EQ(1, m2.rank());
    m2 = fill_square_matrix<double>(10);
    ASSERT_EQ(1, m2.rank());
}

TEST_F(TestMatrix, inverse) {
    Matrix m(2, 2);
    m[0][0] = 3;
    m[0][1] = 1;
    m[1][0] = 1;
    m[1][1] = 3;
    auto m3 = m.inv();
    ASSERT_TRUE(std::abs(0.375 - m3[0][0]) < 0.000001);
    ASSERT_TRUE(std::abs(-0.125 - m3[0][1]) < 0.000001);
    ASSERT_TRUE(std::abs(-0.125 - m3[1][0]) < 0.000001);
    ASSERT_TRUE(std::abs(0.375 - m3[1][1]) < 0.000001);
    Matrix m1(3, 3);
    m1[0][0] = 3;
    m1[0][1] = 1;
    m1[1][0] = 1;
    m1[1][1] = 3;
    m1[2][2] = 3;
    m3 = m1.inv();
    ASSERT_TRUE(std::abs(0.375 - m3[0][0]) < 0.000001);
    ASSERT_TRUE(std::abs(-0.125 - m3[0][1]) < 0.000001);
    ASSERT_TRUE(std::abs(m3[0][2]) < 0.000001);
    ASSERT_TRUE(std::abs(-0.125 - m3[1][0]) < 0.000001);
    ASSERT_TRUE(std::abs(0.375 - m3[1][1]) < 0.000001);
    ASSERT_TRUE(std::abs(m3[1][2]) < 0.000001);
    ASSERT_TRUE(std::abs(m3[2][0]) < 0.000001);
    ASSERT_TRUE(std::abs(m3[2][1]) < 0.000001);
    ASSERT_TRUE(std::abs(1.0 / 3 - m3[2][2]) < 0.000001);
    Matrix m2 = fill_rectangle_matrix<double>(2, 3);
    ASSERT_THROW(m2.inv(), std::runtime_error);
    m2 = fill_rectangle_matrix<double>(3, 2);
    ASSERT_THROW(m2.inv(), std::runtime_error);
    m2 = fill_square_matrix<double>(10);
    ASSERT_THROW(m2.inv(), std::runtime_error);
}

TEST_F(TestMatrix, multiplication) {
    Matrix<int32_t> m1 = fill_rectangle_matrix<int32_t>(3, 2);
    Matrix<int32_t> m2 = fill_rectangle_matrix<int32_t>(2, 3);
    Matrix<int32_t> m3(m1);
    m1 *= m2;
    ASSERT_TRUE(m1 == m3 * m2);
    Matrix<int32_t> result(3, 3);
    result[0][0] = 5;
    result[0][1] = 10;
    result[0][2] = 15;
    result[1][0] = 10;
    result[1][1] = 20;
    result[1][2] = 30;
    result[2][0] = 15;
    result[2][1] = 30;
    result[2][2] = 45;
    ASSERT_TRUE(result == m1);
    ASSERT_THROW(m3 * m1, std::runtime_error);
}

TEST_F(TestMatrix, cast) {
    Matrix<int32_t> m1 = fill_rectangle_matrix<int32_t>(2, 3);
    Vector<int32_t> v = static_cast<Vector<int32_t>>(fill_rectangle_matrix<int32_t>(3, 1));
    auto res = m1 * v;
    ASSERT_EQ(2, res.get_size());
    ASSERT_EQ(14, res[0]);
    ASSERT_EQ(28, res[1]);
    auto m = static_cast<Matrix<int32_t>>(res);
    ASSERT_EQ(14, m[0][0]);
    ASSERT_EQ(28, m[1][0]);
}