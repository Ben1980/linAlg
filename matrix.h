#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <vector>
#include <utility>
#include <cassert>

template<typename T>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");

public:
    Matrix(const std::vector<std::vector<T>> &m) : matrix(m) {
        AssertData(matrix);
    }

    Matrix(std::vector<std::vector<T>> &&m) : matrix(std::move(m)) {
        AssertData(matrix);
    }

    const auto & operator()() const {
        return matrix;
    }

    const auto & operator[](size_t index) const {
        return matrix[index];
    }

    template<typename TT>
    friend inline auto operator*(const Matrix<TT> &lhs, const Matrix<TT> & rhs);

private:
    static void AssertData(const std::vector<std::vector<T>> &m) {
        assert(!m.empty());
        const size_t columnCount = m.front().size();

        std::for_each(m.begin(), m.end(), [columnCount](const auto &element) {
            assert(!element.empty());
            assert(element.size() == columnCount);
        });
    }

    static auto Multiply(const Matrix<T> &lhs, const Matrix<T>& rhs) {
        AssertData(lhs());
        AssertData(rhs());
        assert(lhs().front().size() == rhs().size());

        const size_t lhsRows = lhs().size();
        const size_t rhsColumns = rhs().front().size();
        const size_t lhsColumns = lhs().front().size();

        std::vector<std::vector<T>> C(lhsRows, std::vector<T>(rhsColumns));

        for (size_t i = 0; i < lhsRows; ++i) {
            for (size_t k = 0; k < rhsColumns; ++k) {
                for (size_t j = 0; j < lhsColumns; ++j) {
                    C[i][k] += + lhs[i][j] * rhs[j][k];
                }
            }
        }

        return C;
    }

    const std::vector<std::vector<T>> matrix;
};

template<typename T>
inline auto operator*(const Matrix<T> &lhs, const Matrix<T> & rhs) {
    return Matrix<T>::Multiply(lhs, rhs);
}

namespace TestUtils {
    template<typename T>
    bool ValuesAreEqual(const T &toCheck, const T &expected) {
        static_assert(std::is_arithmetic<T>::value, "T must be numeric");

        if constexpr (std::is_integral_v<T>) {
            return std::fabs(toCheck - expected) == 0;
        }
        else if (std::is_floating_point_v<T>) {
            constexpr T EPSILON = std::numeric_limits<T>::min();
            if (std::fabs(toCheck) < EPSILON && std::fabs(expected) < EPSILON) return true;
            if (std::fabs(toCheck) > EPSILON && std::fabs(expected) < EPSILON) return false;

            return std::fabs(toCheck/expected - 1) <= EPSILON;
        }

        return false;
    }

    template<typename T>
    bool CompareMatrix(const Matrix<T> &toCheck, const Matrix<T> &expected) {
        if(toCheck().empty() || expected().empty()) return false;
        if(toCheck().size() != expected().size()) return false;
        if(toCheck().front().empty() || expected().front().empty()) return false;
        if(toCheck().front().size() != expected().front().size()) return false;

        for(size_t row = 0; row < expected().size(); ++row) {
            for(size_t column = 0; column < expected[row].size(); ++column) {
                const T & data = toCheck[row][column];
                const T & e = expected[row][column];

                if(!ValuesAreEqual(data, e)) {
                    return false;
                }
            }
        }

        return true;
    }
}

TEST_SUITE("Matrix test suite") {
    TEST_CASE ("Matrix Multiplication") {
        const Matrix<double> a = {{{3, 2, 1}, {1, 0, 2}}};
        const Matrix<double> b = {{{1, 2}, {0, 1}, {4, 0}}};

        SUBCASE("c=a*b") {
            //            |1  2|
            //            |0  1|
            //            |4  0|
            //
            // |3  2  1|  |7  8|
            // |1  0  2|  |9  2|

            const Matrix c = a * b;
            const Matrix<double> expected = {{{7, 8}, {9, 2}}};

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=b*a") {
            //         |3  2  1|
            //         |1  0  2|
            //
            // |1  2|  |5  2  5|
            // |0  1|  |1  0  2|
            // |4  0|  |12 8  4|

            const Matrix c = b * a;
            const Matrix<double> expected = {{{5, 2, 5}, {1, 0, 2}, {12, 8, 4}}};

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=a*vec") {
            //            |1|
            //            |0|
            //            |4|
            //
            // |3  2  1|  |7|
            // |1  0  2|  |9|

            const Matrix<double> vec = {{{1}, {0}, {4}}};
            const Matrix c = a * vec;
            const Matrix<double> expected = {{{7}, {9}}};

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=vec*b") {
            //            |1  2|
            //            |0  1|
            //            |4  0|
            //
            // |3  2  1|  |7  8|

            const Matrix<double> vec = {{{3, 2, 1}}};
            const Matrix c = vec * b;
            const Matrix<double> expected = {{{7, 8}}};

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=vec*vec") {
            //            |1|
            //            |0|
            //            |4|
            //
            // |3  2  1|  |7|

            const Matrix<double> vecA = {{{3, 2, 1}}};
            const Matrix<double> vecB = {{{1}, {0}, {4}}};
            const Matrix c = vecA * vecB;

            CHECK(TestUtils::ValuesAreEqual(c().front().front(), 7.));
        }
    }

    TEST_CASE("Matrix Decomposition") {
        //     |1  2  3|   |1  0  0|   |1  2  3|
        // A = |1  1  1| = |1  1  0| * |0 -1 -2|
        //     |3  3  1|   |3  3  1|   |0  0 -2|

        SUBCASE("LU") {
            /*Matrix<double, 3, 3> A = {{{1, 2, 3}, {1, 1, 1}, {3, 3, 1}}};
            LinAlg::Decomposition<double, 3, 3> LU = LinAlg::LUDecomposition(A);

            Matrix<double, 3, 3> expectedL = {{{1, 0, 0}, {1, 1, 0}, {3, 3, 1}}};
            Matrix<double, 3, 3> expectedU = {{{1, 2, 3}, {0, -1, -2}, {0, 0, -2}}};

            CHECK(TestUtils::CompareMatrix(LU.L, expectedL));
            CHECK(TestUtils::CompareMatrix(LU.U, expectedU));*/
            CHECK(false);
        }
    }
}

#endif //LINALG_MATRIX_H
