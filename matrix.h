#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>

template <typename T, size_t nbRows, size_t nbColumns>
class Matrix {
    static_assert(std::is_arithmetic<T>::value, "T must be numeric");

public:
    Matrix(std::vector<std::vector<T>> &&vector) {
        assert(vector.size() == nbRows);
        assert(vector.front().size() == nbColumns);

        auto matrixIter = matrix.begin();
        for(auto && row : vector) {
            std::copy_n(std::make_move_iterator(row.begin()), nbColumns, matrixIter);
            std::advance(matrixIter, nbColumns);
        }
    }

    template<size_t rhsRows, size_t rhsColumns>
    inline Matrix & operator*=(const Matrix<T, rhsRows, rhsColumns> &rhs) {
        assert(nbRows == rhsColumns);
        assert(nbColumns == rhsRows);

        return *this;
    }

private:
    std::array<T, nbRows * nbColumns> matrix;
};

template<typename C, size_t LhsRows, size_t LhsColumns, size_t RhsRows, size_t RhsColumns>
inline auto operator*(Matrix<C, LhsRows, LhsColumns> lhs, const Matrix<C, RhsRows, RhsColumns>& rhs) {
    return lhs *= rhs;
}

namespace TestHelper {
    template <typename T, size_t nbRows, size_t nbColumns>
    bool CompareMatrix(const Matrix<T, nbRows, nbColumns> &toCheck, const Matrix<T, nbRows, nbColumns> &expected) {
        return false;
    }
};

TEST_CASE("Matrix Multiplication") {
    Matrix<double, 2, 3> a = {{{3,2,1},{1,0,2}}};
    Matrix<double, 3, 2> b = {{{1,2},{0,1},{4,0}}};

    Matrix c = a * b;

    Matrix<double, 2, 2> expected = {{{7,8},{9,2}}};

    CHECK(TestHelper::CompareMatrix(expected, expected));
}

#endif //LINALG_MATRIX_H
