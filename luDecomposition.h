#ifndef LINALG_LUDECOMPOSITION_H
#define LINALG_LUDECOMPOSITION_H

#include "matrix.h"

namespace LUDecomposition {
    template<typename T>
    struct Decomposition {
        Matrix<T> L;
        Matrix<T> U;
        Matrix<T> P;
        Matrix<T> Q;
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        std::vector<std::vector<T>> L;
        std::vector<std::vector<T>> U;

        for(size_t column = 0; column < matrix.columns(); ++column) {
            for(size_t row = 0; row < matrix.rows(); ++row) {
                for(size_t col = 0; col < matrix.columns(); ++col) {
                    //L[column, ]
                }
            }
        }

        return {};
    }
}

/*TEST_SUITE("Matrix decomposition test suite") {
    TEST_CASE("Matrix Decomposition") {
        //     |1  2  3|   |1  0  0|   |1  2  3|
        // A = |1  1  1| = |1  1  0| * |0 -1 -2|
        //     |3  3  1|   |3  3  1|   |0  0 -2|

        SUBCASE("LU") {
            Matrix<double> A = {{{1, 2, 3}, {1, 1, 1}, {3, 3, 1}}};
            //LUDecomposition::Decomposition<double> LU = LUDecomposition::Decompose(A);

            Matrix<double> expectedL = {{{1, 0, 0}, {1, 1, 0}, {3, 3, 1}}};
            Matrix<double> expectedU = {{{1, 2, 3}, {0, -1, -2}, {0, 0, -2}}};

            //CHECK(TestUtils::CompareMatrix(LU.L, expectedL));
            //CHECK(TestUtils::CompareMatrix(LU.U, expectedU));
        }
    }
}*/

#endif //LINALG_LUDECOMPOSITION_H
