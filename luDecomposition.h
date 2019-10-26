#ifndef LINALG_LUDECOMPOSITION_H
#define LINALG_LUDECOMPOSITION_H

#include "matrix.h"
#include "matrixFactory.h"

namespace LUDecomposition {
    template<typename T>
    struct Decomposition {
        Matrix<T> L;
        Matrix<T> U;

        Decomposition(const Matrix<T> &matrix) : L(MatrixFactory::IdentityMatrix<T>(matrix.rows())), U(matrix) {}
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        Decomposition<T> decomposition(matrix);

        for(size_t column = 0; column < matrix.columns(); ++column) {
            for(size_t row = column + 1; row < matrix.rows(); ++row) {
                const T & divisor = decomposition.U(column, column);
                if(divisor > 0) {
                    decomposition.L(row, column) = decomposition.U(row, column) / divisor;
                }
                for(size_t col = column; col < matrix.columns(); ++col) {
                    decomposition.U(row, col) -= decomposition.L(row, column) * decomposition.U(column, col);
                }
            }
        }

        return decomposition;
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        static const double EPSILON = 1e-10;
        SUBCASE("LU-Decomposition") {
            //     |1 2 3|        |1 0 0|     |1  2  3|
            // A = |1 1 1| -> L = |1 1 0| U = |0 -1 -2|
            //     |3 3 1|        |3 3 1|     |0  0 -2|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{1, 2, 3, 1, 1, 1, 3, 3, 1}).data()
            };
            LUDecomposition::Decomposition<double> decomposition = LUDecomposition::Decompose(A);

            Matrix test = decomposition.L * decomposition.U;
            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }

        SUBCASE("LU-Decomposition") {
            //     |1 2 3|        |1 0 0|     |1  2  3|
            // A = |1 1 1| -> L = |1 1 0| U = |0 -1 -2|
            //     |3 3 1|        |3 3 1|     |0  0 -2|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}).data()
            };
            LUDecomposition::Decomposition<double> decomposition = LUDecomposition::Decompose(A);

            Matrix test = decomposition.L * decomposition.U;
            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }
    }
}


#endif //LINALG_LUDECOMPOSITION_H
