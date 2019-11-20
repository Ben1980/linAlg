#ifndef LINALG_PIVOTLUDECOMPOSITION_H
#define LINALG_PIVOTLUDECOMPOSITION_H

#include <exception>
#include "matrix.h"
#include "matrixFactory.h"
#include "utils.h"

namespace PivotLUDecomposition {
    template<typename T>
    struct Decomposition {
        Matrix<T> P;
        Matrix<T> L;
        Matrix<T> U;

        Decomposition(const Matrix<T> &matrix) : P(matrix.rows(), matrix.columns()), L(matrix.rows(), matrix.columns()), U(matrix) {}
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        const size_t nbRows = matrix.rows();
        const size_t nbColumns = matrix.columns();
        if(nbRows != nbColumns) {
            throw std::domain_error("Matrix is not square.");
        }

        Decomposition<T> decomposition(matrix);

        decomposition.P = MatrixFactory::IdentityMatrix<T>(nbRows);

        for(size_t k = 0; k < nbRows; ++k) {
            T max = 0;
            size_t pk = 0;
            for(size_t i = k; i < nbRows; ++i) {
                T s = 0;
                for(size_t j = k; j < nbColumns; ++j) {
                    s += std::fabs(decomposition.U(i, j));
                }
                T q = std::fabs(decomposition.U(i, k)) / s;
                if(q > max) {
                    max = q;
                    pk = i;
                }
            }

            if(std::fabs(max) < std::numeric_limits<T>::min()) {
                throw std::domain_error("Pivot has 0 value.");
            }

            if(pk != k) {
                for(size_t j = 0; j < nbColumns; ++j) {
                    std::swap(decomposition.P(k, j), decomposition.P(pk, j));
                    std::swap(decomposition.L(k, j), decomposition.L(pk, j));
                    std::swap(decomposition.U(k, j), decomposition.U(pk, j));
                }
            }

            for(size_t i = k+1; i < nbRows; ++i) {
                decomposition.L(i, k) = decomposition.U(i, k)/decomposition.U(k, k);

                for(size_t j = k; j < nbColumns; ++j) {
                    decomposition.U(i, j) = decomposition.U(i, j) - decomposition.L(i, k) * decomposition.U(k, j);
                }
            }
        }

        for(size_t k = 0; k < nbRows; ++k) {
            decomposition.L(k,k) = 1;
        }

        return decomposition;
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        static const double EPSILON = 1e-10;
        SUBCASE("Pivot-LU-Decomposition Test 1") {
            //     |1 2 3|
            // A = |1 1 1|
            //     |3 3 1|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{1, 2, 3, 1, 1, 1, 3, 3, 1}).data()
            };
            PivotLUDecomposition::Decomposition<double> decomposition = PivotLUDecomposition::Decompose(A);

            Matrix test1 = decomposition.L * decomposition.U;
            
            Matrix<double> P = {
                3, 3, (std::array<double, 9>{0, 0, 1, 1, 0, 0, 0, 1, 0}).data()
            };
            CHECK(TestUtils::CompareMatrix(decomposition.P, P, false, EPSILON));

            Matrix test2 = decomposition.P * A;

            CHECK(TestUtils::CompareMatrix(test1, test2, false, EPSILON));
            
        }

        SUBCASE("Pivot-LU-Decomposition Test 2") {
            //     |2.1  2512 -2516|
            // A = |-1.3  8.8  -7.6|
            //     |0.9   -6.2  4.6|

            Matrix<double> A = {
                3, 3, (std::array<double, 9>{2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}).data()
            };
            PivotLUDecomposition::Decomposition<double> decomposition = PivotLUDecomposition::Decompose(A);

            Matrix test1 = decomposition.L * decomposition.U;

            Matrix<double> P = {
                3, 3, (std::array<double, 9>{0, 0, 1, 1, 0, 0, 0, 1, 0}).data()
            };
            CHECK(TestUtils::CompareMatrix(decomposition.P, P, false, EPSILON));

            Matrix test2 = P*A;

            CHECK(TestUtils::CompareMatrix(test1, test2, EPSILON));
        }
    }
}

#endif //LINALG_PIVOTLUDECOMPOSITION_H