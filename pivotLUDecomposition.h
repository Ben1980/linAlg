#ifndef LINALG_PIVOTLUCOMPOSITION_H
#define LINALG_PIVOTLUCOMPOSITION_H

#include <cassert>
#include "matrix.h"
#include "matrixFactory.h"

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
        assert(nbRows == nbColumns);

        Decomposition<T> decomposition(matrix);

        for(size_t column = 0; column < nbColumns - 1; ++column) {
            T qi = 0;
            size_t index = 0;
            for(size_t row = column; row < nbRows; ++row) {
                T si = 0;
                for(size_t col = column; col < nbColumns; ++col) {
                    si += std::fabs(decomposition.U(row, col));
                }

                const T tmpQi = std::fabs(decomposition.U(row, column)) / si;
                if(tmpQi > qi) {
                    qi = tmpQi;
                    index = row;
                }
            }

            if(index != column) {
                decomposition.P(index, column) = 1;
                for(size_t col = 0; col < nbColumns; ++col) {
                    std::swap(decomposition.U(column, col), decomposition.U(index, col));
                    std::swap(decomposition.L(column, col), decomposition.L(index, col));
                }
            }

            for(size_t row = column + 1; row < nbRows; ++row) {
                const T & divisor = decomposition.U(column, column);
                assert(std::fabs(divisor) > std::numeric_limits<T>::min()); //a_ii != 0 is necessary because of pivoting with diaognal strategy
                
                decomposition.L(row, column) = decomposition.U(row, column) / divisor;

                for(size_t col = column; col < nbColumns; ++col) {
                    decomposition.U(row, col) -= decomposition.L(row, column) * decomposition.U(column, col);
                }
            }
        }

        for(size_t column = 0; column < nbColumns; ++column) {
            decomposition.L(column, column) = 1;
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
            Matrix test2 = decomposition.P * A;
            CHECK(TestUtils::CompareMatrix(test1, test2, EPSILON));
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
            Matrix test2 = decomposition.P * A;
            CHECK(TestUtils::CompareMatrix(test1, test2, EPSILON));
        }
    }
}

#endif //LINALG_PIVOTCOMPOSITION_H