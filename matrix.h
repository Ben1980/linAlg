#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <memory>
#include <exception>

template<typename T>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");

public:
    ~Matrix() = default;
    Matrix(size_t rows, size_t columns, T *m) : nbRows(rows), nbColumns(columns) {
        const size_t size = nbRows*nbColumns;

        matrix = std::make_unique<T[]>(size);
        std::copy(m, m + size, matrix.get());

        AssertData(*this);
    }
    Matrix(size_t rows, size_t columns) : nbRows(rows), nbColumns(columns) {
        const int size = nbRows * nbColumns;
        matrix = std::make_unique<T[]>(size);

        std::fill(matrix.get(), matrix.get() + size, 0);
        AssertData(*this);
    }
    Matrix(const Matrix<T> &m) : nbRows(m.nbRows), nbColumns(m.nbColumns) {
        const int size = nbRows * nbColumns;

        matrix = std::make_unique<T[]>(size);
        std::copy(m.matrix.get(), m.matrix.get() + size, matrix.get());
    }
    Matrix(Matrix<T> &&m) : nbRows(std::move(m.nbRows)), nbColumns(std::move(m.nbColumns)) {
        matrix.swap(m.matrix);

        m.nbRows = 0;
        m.nbColumns = 0;
        m.matrix.release();
    }
    Matrix<T> & operator=(const Matrix<T> &m){
        Matrix tmp(m);

        nbRows = tmp.nbRows;
        nbColumns = tmp.nbColumns;
        matrix.reset(tmp.matrix.get());

        return *this;
    }
    Matrix<T> & operator=(Matrix<T> &&m){
        Matrix tmp(std::move(m));

        std::swap(tmp.nbRows, nbRows);
        std::swap(tmp.nbColumns, nbColumns);
        matrix.swap(tmp.matrix);

        return *this;
    }
    
    const T & operator()(size_t row, size_t column) const {
        return matrix[row*nbColumns + column];
    }

    T & operator()(size_t row, size_t column) {
        return matrix[row*nbColumns + column];
    }

    [[nodiscard]] size_t rows() const {
        return nbRows;
    }

    [[nodiscard]] size_t columns() const {
        return nbColumns;
    }

    template<typename U>
    friend Matrix<U> operator*(const Matrix<U> &lhs, const Matrix<U> & rhs);

private:
    static void AssertData(const Matrix<T> &m) {
        if(m.nbRows < 1 || m.nbColumns < 1) {
          throw std::domain_error("Invalid defined matrix.");
        }
        if(m.nbRows != m.nbColumns) {
          throw std::domain_error("Matrix is not square.");
        }
    }

    size_t nbRows{0};
    size_t nbColumns{0};

    std::unique_ptr<T[]> matrix{nullptr};
};

template<typename U>
Matrix<U> operator*(const Matrix<U> &lhs, const Matrix<U> & rhs) {
    Matrix<U>::AssertData(lhs);
    Matrix<U>::AssertData(rhs);
    if(lhs.rows() != rhs.rows()) {
      throw std::domain_error("Matrices have unequal size.");
    }

    const size_t lhsRows = lhs.rows();
    const size_t rhsColumns = rhs.columns();
    const size_t lhsColumns = lhs.columns();

    Matrix<U> C(lhsRows, rhsColumns);
    
    for (size_t i = 0; i < lhsRows; ++i) {
        for (size_t k = 0; k < rhsColumns; ++k) {
            for (size_t j = 0; j < lhsColumns; ++j) {
                C(i, k) += lhs(i, j) * rhs(j, k);
            }
        }
    }

    return C;
}

namespace TestUtils {
    template<typename T>
    bool ValuesAreEqual(const T &toCheck, const T &expected, T epsilon = std::numeric_limits<T>::min()) {
        static_assert(std::is_arithmetic<T>::value, "T must be numeric");

        if constexpr (std::is_integral_v<T>) {
            return std::fabs(toCheck - expected) == 0;
        }
        else if (std::is_floating_point_v<T>) {
            if (std::fabs(toCheck) < epsilon && std::fabs(expected) < epsilon) return true;
            if (std::fabs(toCheck) > epsilon && std::fabs(expected) < epsilon) return false;

            return std::fabs(toCheck/expected - 1) <= epsilon;
        }

        return false;
    }

    template<typename T>
    bool CompareMatrix(const Matrix<T> &toCheck, const Matrix<T> &expected, T epsilon = std::numeric_limits<T>::min()) {
        if(toCheck.rows() == 0 || expected.rows() == 0) return false;
        if(toCheck.columns() == 0 || expected.columns() == 0) return false;
        if(toCheck.rows() != expected.rows()) return false;
        if(toCheck.columns() != expected.columns()) return false;

        bool equalData = true;
        for(size_t row = 0; row < expected.rows(); ++row) {
            for(size_t column = 0; column < expected.columns(); ++column) {
                const T & data = toCheck(row, column);
                const T & e = expected(row, column);

                equalData = equalData && ValuesAreEqual(data, e, epsilon);
            }
        }

        if(!equalData) {
            std::vector<T> check(toCheck.rows()*toCheck.columns()), expect(expected.rows()*expected.columns());

            for(size_t row = 0; row < expected.rows(); ++row) {
                for(size_t column = 0; column < expected.columns(); ++column) {
                    check[row*toCheck.columns() + column] = toCheck(row, column);
                    expect[row*expected.columns() + column] = expected(row, column);
                }
            }

            fmt::print("To Check = {}\n", fmt::join(check, " "));
            fmt::print("Expected = {}\n", fmt::join(expect, " "));
        }
        
        return equalData;
    }
}

TEST_SUITE("Matrix test suite") {
    TEST_CASE ("Matrix Multiplication") {
        const Matrix<double> a = {
            2, 3, (std::array<double, 6>{3, 2, 1, 1, 0, 2}).data()
        };
        const Matrix<double> b = {
            3, 2, (std::array<double, 6>{1, 2, 0, 1, 4, 0}).data()
        };

        SUBCASE("c=a*b") {
            //            |1  2|
            //            |0  1|
            //            |4  0|
            //
            // |3  2  1|  |7  8|
            // |1  0  2|  |9  2|

            Matrix c = a * b;
            const Matrix<double> expected = {
                2, 2, (std::array<double, 4>{7, 8, 9, 2}).data()
            };

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=b*a") {
            //         |3  2  1|
            //         |1  0  2|
            //
            // |1  2|  |5  2  5|
            // |0  1|  |1  0  2|
            // |4  0|  |12 8  4|

            Matrix c = b * a;
            const Matrix<double> expected = {
                3, 3, (std::array<double, 9>{5, 2, 5, 1, 0, 2, 12, 8, 4}).data()
            };

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=a*vec") {
            //            |1|
            //            |0|
            //            |4|
            //
            // |3  2  1|  |7|
            // |1  0  2|  |9|

            const Matrix<double> vec = {
                3, 1, (std::array<double, 3>{1, 0, 4}).data()
            };
            Matrix c = a * vec;
            const Matrix<double> expected = {
                2, 1, (std::array<double, 2>{7, 9}).data()
            };

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=vec*b") {
            //            |1  2|
            //            |0  1|
            //            |4  0|
            //
            // |3  2  1|  |7  8|

            const Matrix<double> vec = {
                1, 3, (std::array<double, 3>{3, 2, 1}).data()
            };
            Matrix c = vec * b;
            const Matrix<double> expected = {
                1, 2, (std::array<double, 2>{7, 8}).data()
            };

            CHECK(TestUtils::CompareMatrix(c, expected));
        }

        SUBCASE("c=vec*vec") {
            //            |1|
            //            |0|
            //            |4|
            //
            // |3  2  1|  |7|

            const Matrix<double> vecA = {
                1, 3, (std::array<double, 3>{3, 2, 1}).data()
            };
            const Matrix<double> vecB = {
                3, 1, (std::array<double, 3>{1, 0, 4}).data()
            };
            Matrix c = vecA * vecB;

            CHECK(TestUtils::ValuesAreEqual(c(0, 0), 7.));
        }
    }
}

#endif //LINALG_MATRIX_H
