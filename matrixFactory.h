#ifndef LINALG_MATRIXFACTORY_H
#define LINALG_MATRIXFACTORY_H

#include <exception>
#include <random>
#include <vector>
#include <array>
#include "matrix.h"

namespace MatrixFactory {
    template <typename T>
    Matrix<T> IdentityMatrix(size_t size) {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");
        if(size == 0) {
            throw std::domain_error("Invalid defined matrix.");
        }

        const size_t arrSize = size*size;
        std::vector<T> vec(arrSize);
        for(size_t index = 0; index < size; ++index) {
            vec[index*size + index] = 1;
        }

        Matrix<T> matrix(size, size, &vec[0]);
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
    Matrix<T> RandomMatrix(const Range<T> &range) {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");
        if(range.rows == 0 || range.columns == 0) {
            throw std::domain_error("Invalid defined matrix.");
        }

        Matrix<T> matrix(range.rows, range.columns);

        std::mt19937 mt(std::random_device{}());
        std::uniform_real_distribution<> real_dist(range.from, range.to);
        const auto gen = std::bind(std::ref(real_dist), std::ref(mt));

        for(size_t row = 0; row < range.rows; ++row) {
            for(size_t column = 0; column < range.columns; ++column) {
                matrix(row, column) = gen();
            }
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
        Matrix<int> expected = {
            4, 4, (std::array<int, 16>{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}).data()
        };

        CHECK(TestUtils::CompareMatrix(identity, expected));
    }

    TEST_CASE ("Random Matrix") {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({4, 4, 1, 5});

        for(size_t row = 0; row < A.rows(); ++row) {
            for(size_t column = 0; column < A.columns(); ++column) {
                const bool isValid = A(row, column) >= 1 && A(row, column) <= 5;
                CHECK(isValid);
            }
        }
    }
}

#endif //LINALG_MATRIXFACTORY_H
