#ifndef LINALG_CHOLESKYDECOMPOSITION_H
#define LINALG_CHOLESKYDECOMPOSITION_H

#include <exception>
#include "matrix.h"
#include "matrixFactory.h"

namespace CholeskyDecomposition {
    template<typename T>
    Matrix<T> SimplifySymmetricMatrix(Matrix<T> matrix) {
        const size_t nbRows = matrix.rows();
        const size_t nbColumns = matrix.columns();

        for(size_t row = 0; row < nbRows; ++row) {
            for(size_t column = row + 1; column < nbColumns; ++column) {
                matrix(row, column) = 0;
            }
        }

        return matrix;
    }

    template<typename T>
    Matrix<T> Decompose(const Matrix<T> &matrix) {
        const size_t nbRows = matrix.rows();
        const size_t nbColumns = matrix.columns();
        if(nbRows != nbColumns) {
            throw std::domain_error("Matrix is not square.");
        }

        Matrix<T> L = SimplifySymmetricMatrix(matrix);
        
        for(size_t k = 0; k < nbColumns; ++k) {
            const T & a_kk = L(k, k);
            
            if(a_kk > 0) {
                L(k, k) = std::sqrt(a_kk);

                for(size_t i = k + 1; i < nbRows; ++i) {
                    L(i, k) /= L(k, k);
                    
                    for(size_t j = k + 1; j <= i; ++j) {
                        L(i, j) -= L(i, k) * L(j, k);
                    }
                }
            }
            else {
                throw std::domain_error("Matrix is not positive definit.");
            }
        }
        
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

            CHECK(TestUtils::CompareMatrix(test, A, false, EPSILON));
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
            CHECK(TestUtils::CompareMatrix(test, A, false, EPSILON));
        }

        SUBCASE("Cholesky-Decomposition Test 3, Error") {
            //     |0  7  3|
            // A = |7 11  2|
            //     |3  2  6|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{0, 7, 3, 7, 11, 2, 3, 2, 6}).data()
            };
            CHECK_THROWS_AS(CholeskyDecomposition::Decompose(A), std::domain_error);
        }
    }
}

#endif //LINALG_CHOLESKYDECOMPOSITION_H