#ifndef LINALG_CHOLESKYDECOMPOSITION_H
#define LINALG_CHOLESKYDECOMPOSITION_H

#include <exception>
#include "matrix.h"
#include "matrixFactory.h"

namespace CholeskyDecomposition {
    template<typename T>
    Matrix<T> Decompose(const Matrix<T> &matrix) {
        const size_t nbRows = matrix.rows();
        const size_t nbColumns = matrix.columns();
        if(nbRows != nbColumns) {
            throw std::domain_error("Matrix is not square.");
        }

        Matrix<T> L(nbRows, nbColumns);
        for(size_t i = 0; i < nbRows; i++) {
            for(size_t j = 0; j < nbColumns; j++) {
                
            }
        }

        /*for (size_t k = 0; k < nbColumns; ++k) {
            const T & a_kk = L(k, k);

            if(a_kk < std::numeric_limits<T>::min()) {
                throw std::domain_error("Matrix is not positive definite.");
            }

            const T l_kk = std::sqrt(a_kk);
            for (size_t i = k + 1; i < nbRows; ++i) {
                const T l_ik = L(i, k) / l_kk;
                
                for (size_t j = k + 1; j < nbColumns; ++j) {
                    L(i, j) = L(i, j) - l_ik*L(j, k) / l_kk;
                }
            }
        }*/

        return L;
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        static const double EPSILON = 1e-10;
        SUBCASE("Cholesky-Decomposition Test 1") {
            //     |5  7  3|
            // A = |7 11  2|
            //     |3  2  6|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{5, 7, 3, 7, 11, 2, 3, 2, 6}).data()
            };
            Matrix L = CholeskyDecomposition::Decompose(A);

            Matrix test = L * L.transpose();

            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }

        SUBCASE("Cholesky-Decomposition Test 2") {
            //     | 9   3   -6   12|
            // A = | 3   26  -7  -11|
            //     |-6  -7    9    7|
            //     |12  -11   7   65|

            Matrix<double> A = {
                4, 4, (std::array<double, 16>{9, 3, -6, 12, 3, 26, -7, -11, -6, -7, 9, 7, 12, -11, 7, 65}).data()
            };
            Matrix L = CholeskyDecomposition::Decompose(A);

            Matrix test = L * L.transpose();
            CHECK(TestUtils::CompareMatrix(test, A, EPSILON));
        }
    }
}

#endif //LINALG_CHOLESKYDECOMPOSITION_H