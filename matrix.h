#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <cmath>
#include <vector>
#include <type_traits>
#include <cassert>
#include <algorithm>

template<typename T>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");

public:
#ifdef _ASMATRIX
    Matrix(size_t rows, size_t columns, T *m) : nbRows(rows), nbColumns(columns) {
        matrix.resize(nbRows);
        for(size_t rowIndex = 0; rowIndex < nbRows; ++rowIndex) {
            auto & row = matrix[rowIndex];
            row.reserve(nbColumns);
            for(size_t column = 0; column < nbColumns; ++column) {
                row.push_back(m[rowIndex*nbColumns + column]);
            }
        }
        AssertData(*this);
    }
    Matrix(size_t rows, size_t columns) : nbRows(rows), nbColumns(columns), matrix(std::vector<std::vector<T>>(rows, std::vector<T>(columns))) {
        AssertData(*this);
    }
#elif _ASARRAY
    Matrix() : matrix(nullptr) {}
    Matrix(size_t nbRows, size_t nbColumns, T *m) : nbRows(nbRows), nbColumns(nbColumns) {
        const size_t size = nbRows*nbColumns;

        matrix = new T[size];
        std::copy(m, m + size, matrix);

        AssertData(*this);
    }
    Matrix(size_t nbRows, size_t nbColumns) : nbRows(nbRows), nbColumns(nbColumns) {
        const int size = nbRows * nbColumns;
        matrix = new T[size];
        for(size_t index = 0; index < size; ++index) {
            matrix[index] = 0;
        }
        AssertData(*this);
    }
    ~Matrix() {
        delete [] matrix;
        matrix = nullptr;
    }
    Matrix(const Matrix &rhs) {
        if(this != &rhs) {
            nbRows = rhs.nbRows;
            nbColumns = rhs.nbColumns;

            if(matrix) {
                delete [] matrix;
            }

            const size_t size = nbRows*nbColumns;
            matrix = new T[size];
            for(size_t index = 0; index < size; ++index) {
                matrix[index] = rhs.matrix[index];
            }
        }
    }
    Matrix(Matrix &&rhs) {
        if(this != &rhs) {
            nbRows = std::move(rhs.nbRows);
            nbColumns = std::move(rhs.nbColumns);
            std::swap(matrix, rhs.matrix);
        }
    }
    Matrix & operator=(const Matrix &rhs) {
        if(this != &rhs) {
            Matrix tmp(rhs);
            swap(tmp, *this);
        }
        return *this;
    }
    Matrix & operator=(Matrix &&rhs) {
        if(this != &rhs) {
            Matrix tmp(rhs);
            swap(tmp, *this);
        }
        return *this;
    }
#else
    Matrix(size_t nbRows, size_t nbColumns, T *m) : nbRows(nbRows), nbColumns(nbColumns) {
        const int size = nbRows * nbColumns;
        matrix.reserve(size);
        for(size_t index = 0; index < size; ++index) {
            matrix.push_back(m[index]);
        }
        AssertData(*this);
    }
    Matrix(size_t nbRows, size_t nbColumns) : nbRows(nbRows), nbColumns(nbColumns), matrix(std::vector<T>(nbRows*nbColumns)) {
        AssertData(*this);
    }
#endif
    
    const T & operator()(size_t row, size_t column) const {
#ifdef _ASMATRIX       
        return matrix[row][column];
#else
        return matrix[row*nbColumns + column];
#endif
    }

    T & operator()(size_t row, size_t column) {
#ifdef _ASMATRIX       
        return matrix[row][column];
#else
        return matrix[row*nbColumns + column];
#endif
    }

    [[nodiscard]] size_t rows() const {
        return nbRows;
    }

    [[nodiscard]] size_t columns() const {
        return nbColumns;
    }

    template<typename U>
    friend inline auto operator*(const Matrix<U> &lhs, const Matrix<U> & rhs);

private:
    static void AssertData(const Matrix<T> &m) {
#ifdef _ASMATRIX 
        assert(!m.matrix.empty());
        assert(!m.matrix.front().empty());

        for(const auto & row : m.matrix) {
            assert(row.size() == m.nbColumns);
        }
#elif _ASARRAY
        assert(m.nbRows > 0);
        assert(m.nbColumns > 0);
        assert(m.matrix != nullptr);
#else
        assert(!m.matrix.empty());
        assert(m.matrix.size() == m.nbRows*m.nbColumns);
#endif
    }

    static Matrix<T> Multiply(const Matrix<T> &lhs, const Matrix<T>& rhs) {
        AssertData(lhs);
        AssertData(rhs);
        assert(lhs.columns() == rhs.rows());

        const size_t lhsRows = lhs.rows();
        const size_t rhsColumns = rhs.columns();
        const size_t lhsColumns = lhs.columns();

        Matrix<T> C(lhsRows, rhsColumns);
        
        for (size_t i = 0; i < lhsRows; ++i) {
            for (size_t k = 0; k < rhsColumns; ++k) {
                for (size_t j = 0; j < lhsColumns; ++j) {
                    C(i, k) += lhs(i, j) * rhs(j, k);
                }
            }
        }

        return C;
    }

    static void swap(Matrix &lhs, Matrix &rhs) {
        std::swap(lhs.nbRows, rhs.nbRows);
        std::swap(lhs.nbColumns, rhs.nbColumns);
        std::swap(lhs.matrix, rhs.matrix);
    }

    size_t nbRows{0};
    size_t nbColumns{0};

#ifdef _ASMATRIX
    std::vector<std::vector<T>> matrix;
#elif _ASARRAY
    T * matrix{nullptr};
#else
    std::vector<T> matrix;
#endif
};

template<typename U>
inline auto operator*(const Matrix<U> &lhs, const Matrix<U> & rhs) {
    return Matrix<U>::Multiply(lhs, rhs);
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
        if(toCheck.rows() == 0 || expected.rows() == 0) return false;
        if(toCheck.columns() == 0 || expected.columns() == 0) return false;
        if(toCheck.rows() != expected.rows()) return false;
        if(toCheck.columns() != expected.columns()) return false;

        bool equalData = true;
        for(size_t row = 0; row < expected.rows(); ++row) {
            for(size_t column = 0; column < expected.columns(); ++column) {
                const T & data = toCheck(row, column);
                const T & e = expected(row, column);

                equalData = equalData && ValuesAreEqual(data, e);
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
