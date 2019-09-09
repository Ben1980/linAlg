#ifndef LINALG_GAUSS_H
#define LINALG_GAUSS_H

#include "matrix.h"

namespace Gauss {
    template<typename T>
    Matrix<T> Solve(const Matrix<T> &matrix, const std::vector<T> &r) {
        assert(matrix.rows() == matrix.columns());
        assert(r.size() == matrix.rows());



        return {};
    }
}

TEST_SUITE("Matrix solve test suite") {
    TEST_CASE("Matrix Decomposition") {
        //     |2  3 -5|-10|   |2  3 -5|-10|
        // A = |4  8 -3|-19| = |0  2  7| 1 |
        //     |-6 1  4|-11|   |0  0-46|-46|

        SUBCASE("Gauss") {
            Matrix<double> A = {{{2, 3, -5}, {4, 8, -3}, {-6, 1, 4}}};
            Matrix<double> B = Gauss::Solve(A, {-10, -19, -11});

            Matrix<double> expected = {{{2, 3, -5}, {0, 2, 7}, {0, 0, -46}}};

            CHECK(TestUtils::CompareMatrix(B, expected));
        }
    }
}

#endif //LINALG_GAUSS_H
