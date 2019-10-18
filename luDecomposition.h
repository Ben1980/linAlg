#ifndef LINALG_LUDECOMPOSITION_H
#define LINALG_LUDECOMPOSITION_H

#include "matrix.h"

namespace LUDecomposition {
    template<typename T>
    struct Decomposition {
        Matrix<T> L;
        Matrix<T> U;

        Decomposition(const Matrix<T> &matrix) : L(matrix.rows(), matrix.columns()), U(matrix) {}
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        Decomposition<T> decomposition(matrix);

        for(size_t column = 0; column < matrix.columns(); ++column) {
            for(size_t row = column+1; row < matrix.rows(); ++row) {
                const T & val = decomposition.U(column, column);
                if(val > 0) {
                    decomposition.L(row, column) = decomposition.U(row, column) / decomposition.U(column, column);
                }
                for(size_t col = column; col < matrix.columns(); ++col) {
                    decomposition.U(row, col) = decomposition.U(row, col) - decomposition.L(row, column) * decomposition.U(column, col);
                }
            }
        }

        return decomposition;
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        //     |2.1  2512 -2512|        |2.1  2512     -2512|     |0        0       0|
        // A = |-1.3  8.8  -7.6| -> L = |0    1563.9 -1565.1| R = |-0.61905 0       0|
        //     |0.9  -6.2   4.6|        |0    0         -0.7|     |0.42857 -0.69237 0|

        SUBCASE("LU-Decomposition") {
            Matrix<double> A = {
                3, 3, (std::array<double, 9>{2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}).data()
            };
            LUDecomposition::Decomposition<double> decomposition = LUDecomposition::Decompose(A);

            Matrix<double> expectedU = {
                3, 3, (std::array<double, 9>{2.1, 2512, -2512, 0, 1563.9, -1565.1, 0, 0, -0.7}).data()
            };

            Matrix<double> expectedL = {
                3, 3, (std::array<double, 9>{0, 0, 0, -0.61905, 0, 0, 0.42857, -0.69237, 0}).data()
            };

            CHECK(TestUtils::CompareMatrix(decomposition.L, expectedL));
            CHECK(TestUtils::CompareMatrix(decomposition.U, expectedU));
        }
    }
}


#endif //LINALG_LUDECOMPOSITION_H
