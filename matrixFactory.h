#ifndef LINALG_MATRIXFACTORY_H
#define LINALG_MATRIXFACTORY_H

#include <random>
#include <cassert>
#include "matrix.h"

namespace MatrixFactory {
    template <typename T>
    auto IdentityMatrix(size_t size) {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");
        assert(size > 0);

        std::vector<std::vector<T>> matrix(size, std::vector<T>(size));

        for(size_t index = 0; index < size; ++index) {
            matrix[index][index] = 1;
        }

        return matrix;
    }

    template<typename T>
    struct Range {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");

        size_t rows{0};
        size_t columns{0};
        T from;
        T to;
    };

    template <typename T>
    auto RandomMatrix(const Range<T> &range) {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");
        assert(range.rows > 0);
        assert(range.columns > 0);

        std::vector<std::vector<T>> matrix(range.rows, std::vector<T>(range.columns));

        std::mt19937 mt(std::random_device{}());
        std::uniform_real_distribution<> real_dist(range.from, range.to);
        const auto gen = std::bind(std::ref(real_dist), std::ref(mt));

        for(auto & row : matrix) {
            std::generate(row.begin(), row.end(), [gen]() -> auto {
                return gen();
            });
        }

        return matrix;
    }
}

TEST_SUITE("MatrixFactory test suite") {
    TEST_CASE ("Identity Matrix") {
        // |1  0  0  0|
        // |0  1  0  0|
        // |0  0  1  0|
        // |0  0  0  1|

        Matrix<int> identity = MatrixFactory::IdentityMatrix<int>(4);
        Matrix<int> expected = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}};

        CHECK(TestUtils::CompareMatrix(identity, expected));
    }

    TEST_CASE ("Random Matrix") {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({4, 4, 1, 5});

        for(const auto &row : A()) {
            for(const auto &element : row) {
                const bool isValid = element >= 1 && element <= 5;
                CHECK(isValid);
            }
        }
    }
}

#endif //LINALG_MATRIXFACTORY_H
