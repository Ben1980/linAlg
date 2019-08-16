#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <fmt/format.h>
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>

template <typename T, size_t rows, size_t columns>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");

public:
    constexpr Matrix(std::array<T, rows * columns> &&array) : matrix(array) { }

    constexpr Matrix(std::vector<std::vector<T>> &&vector) {
        assert(vector.size() == rows);
        assert(vector.front().size() == columns);

        auto matrixIter = matrix.begin();
        for(auto && row : vector) {
            std::copy_n(std::make_move_iterator(row.begin()), columns, matrixIter);
            std::advance(matrixIter, columns);
        }
    }

    [[nodiscard]] const std::array<T, rows * columns> & getRawData() const { return matrix; }

private:
    std::array<T, rows * columns> matrix;
};

template<typename T, size_t lhsRows, size_t lhsColumns, size_t rhsRows, size_t rhsColumns>
constexpr inline auto operator*(Matrix<T, lhsRows, lhsColumns> lhs, const Matrix<T, rhsRows, rhsColumns>& rhs) {
    assert(lhsRows == rhsColumns);
    assert(lhsColumns == rhsRows);

    const auto &A = lhs.getRawData();
    const auto &B = rhs.getRawData();

    constexpr size_t size = lhsRows * rhsColumns;
    std::array<T, size> C{};

    for (size_t i = 0; i < lhsRows; ++i) {
        for (size_t j = 0; j < rhsColumns; ++j) {
            T value = 0;
            for (size_t k = 0; k < lhsColumns; ++k) {
                const T a = A[k + i * lhsColumns], b = B[k * rhsColumns + j];
                value += a * b;
            }

            C[j + i * rhsColumns] = value;
        }
    }

    return Matrix<T, lhsRows, rhsColumns>(std::move(C));
}

namespace Helper {
    template<typename T, size_t rows, size_t columns>
    constexpr void PrintMatrix(const Matrix<T, rows, columns> &m) {
        const auto & rawData = m.getRawData();
        auto begin = rawData.begin();
        auto end = rawData.begin() + columns;

        for(size_t row = 0; row < rows; ++row) {
            const std::vector<T> rowVec(begin, end);
            fmt::print("| {:^5} |\n", fmt::join(rowVec, ""));

            std::advance(begin, columns);
            std::advance(end, columns);
        }
    }
}

namespace TestHelper {
    template<typename T>
    constexpr bool valuesAreEqual(const T &toCheck, const T &expected) {
        static_assert(std::is_arithmetic<T>::value, "T must be numeric");
        constexpr T EPSILON = std::numeric_limits<T>::min();

        if(std::fabs(toCheck) < EPSILON && std::fabs(expected) < EPSILON) return true;
        if(std::fabs(toCheck) > EPSILON && std::fabs(expected) < EPSILON) return false;

        return std::fabs(toCheck/expected - 1) <= EPSILON;
    }

    template<typename T, size_t lhsRows, size_t lhsColumns, size_t rhsRows, size_t rhsColumns>
    constexpr bool CompareMatrix(const Matrix<T, lhsRows, lhsColumns> &toCheck, const Matrix<T, rhsRows, rhsColumns> &expected) {
        static_assert(std::is_arithmetic<T>::value, "T must be numeric");

        if(lhsRows != rhsRows || lhsColumns != rhsColumns) return false;

        const auto & checkRawData = toCheck.getRawData();
        const auto & expectedRawData = expected.getRawData();

        for(size_t index = 0; index < expectedRawData.size(); ++index) {
            if(!valuesAreEqual(checkRawData[index], expectedRawData[index])) {
                return false;
            }
        }

        return true;
    }
}

TEST_SUITE("Matrix operations test suite") {
    TEST_CASE ("Matrix Multiplication, c=a*b") {
        //            |1  2|
        //            |0  1|
        //            |4  0|
        //
        // |3  2  1|  |7  8|
        // |1  0  2|  |9  2|

        Matrix<double, 2, 3> a = {{{3, 2, 1}, {1, 0, 2}}};
        Matrix<double, 3, 2> b = {{{1, 2}, {0, 1}, {4, 0}}};

        Matrix c = a * b;

        Matrix<double, 2, 2> expected = {{{7, 8}, {9, 2}}};

        CHECK(TestHelper::CompareMatrix(c, expected));
    }

    TEST_CASE ("Matrix Multiplication, c=b*a") {
        //         |3  2  1|
        //         |1  0  2|
        //
        // |1  2|  |5  2  5|
        // |0  1|  |1  0  2|
        // |4  0|  |12 8  4|

        Matrix<double, 2, 3> a = {{{3, 2, 1}, {1, 0, 2}}};
        Matrix<double, 3, 2> b = {{{1, 2}, {0, 1}, {4, 0}}};

        Matrix c = b * a;

        Matrix<double, 3, 3> expected = {{{5, 2, 5}, {1, 0, 2}, {12, 8, 4}}};

        CHECK(TestHelper::CompareMatrix(c, expected));
    }
}

#endif //LINALG_MATRIX_H
