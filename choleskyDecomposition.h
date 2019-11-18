#ifndef LINALG_CHOLESKYDECOMPOSITION_H
#define LINALG_CHOLESKYDECOMPOSITION_H

#include <exception>
#include "matrix.h"
#include "matrixFactory.h"

namespace CholeskyDecomposition {
    template<typename T>
    struct Decomposition {
        Matrix<T> L;
        Matrix<T> U;

        Decomposition(const Matrix<T> &matrix) : L(MatrixFactory::IdentityMatrix<T>(matrix.rows())), U(matrix) {}
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        throw std::invalid_argument("Not implemented!");
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        static const double EPSILON = 1e-10;
        SUBCASE("Cholesky-Decomposition Test 1") {
            //     |1 2 3|
            // A = |1 1 1|
            //     |3 3 1|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{1, 2, 3, 1, 1, 1, 3, 3, 1}).data()
            };
            CholeskyDecomposition::Decomposition<double> decomposition = CholeskyDecomposition::Decompose(A);

            Matrix test = decomposition.L * decomposition.U;

            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }

        SUBCASE("Cholesky-Decomposition Test 2") {
            //     |2.1  2512 -2516|
            // A = |-1.3  8.8  -7.6|
            //     |0.9   -6.2  4.6|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}).data()
            };
            CholeskyDecomposition::Decomposition<double> decomposition = CholeskyDecomposition::Decompose(A);

            Matrix test = decomposition.L * decomposition.U;
            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }
    }
}

#endif //LINALG_CHOLESKYDECOMPOSITION_H