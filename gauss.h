#ifndef LINALG_GAUSS_H
#define LINALG_GAUSS_H

#include "matrix.h"

namespace Gauss {
    template<typename T>
    struct Decomposition {
        Matrix<T> L;
        Matrix<T> U;
    };

    template<typename T>
    Decomposition<T> Decompose(const Matrix<T> &matrix) {
        assert(matrix.rows() == matrix.columns());



        return {};
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        //     |2.1  2512 -2512|        |2.1  2512     -2512|     |0        0       0|
        // A = |-1.3  8.8  -7.6| -> L = |0    1563.9 -1565.1| R = |-0.61905 0       0|
        //     |0.9  -6.2   4.6|        |0    0         -0.7|     |0.42857 -0.69237 0|

        SUBCASE("Gauss") {
            Matrix<double> A = {
#ifdef _ASMATRIX
                {{2.1, 2512, -2516}, {-1.3, 8.8, -7.6}, {0.9, -6.2, 4.6}}
#elif _ASARRAY
                3, 3, (std::array<double, 9>{2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}).data()
#else
                3, 3, {2.1, 2512, -2516, -1.3, 8.8, -7.6, 0.9, -6.2, 4.6}
#endif
            };
            Gauss::Decomposition<double> decomposition = Gauss::Decompose(A);

            Matrix<double> expectedU = {
#ifdef _ASMATRIX
                {{2.1, 2512, -2512}, {0, 1563.9, -1565.1}, {0, 0, -0.7}}
#elif _ASARRAY
                3, 3, (std::array<double, 9>{2.1, 2512, -2512, 0, 1563.9, -1565.1, 0, 0, -0.7}).data()
#else
                3, 3, {2.1, 2512, -2512, 0, 1563.9, -1565.1, 0, 0, -0.7}
#endif
            };

            Matrix<double> expectedL = {
#ifdef _ASMATRIX
                {{0, 0, 0}, {-0.61905, 0, 0}, {0.42857, -0.69237, 0}}
#elif _ASARRAY
                3, 3, (std::array<double, 9>{0, 0, 0, -0.61905, 0, 0, 0.42857, -0.69237, 0}).data()
#else
                3, 3, {0, 0, 0, -0.61905, 0, 0, 0.42857, -0.69237, 0}
#endif
            };

            CHECK(TestUtils::CompareMatrix(decomposition.L, expectedL));
            CHECK(TestUtils::CompareMatrix(decomposition.U, expectedU));
        }
    }
}

#endif //LINALG_GAUSS_H
