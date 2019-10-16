#ifndef LINALG_MATRIXFACTORY_H
#define LINALG_MATRIXFACTORY_H

#include <random>
#include <cassert>
#include <functional>

namespace MatrixFactory {
    template <typename T>
    auto IdentityMatrix(size_t size) {
        static_assert(std::is_arithmetic<T>::value, "C must be numeric");
        assert(size > 0);

#ifdef _ASMATRIX
        std::vector<std::vector<T>> matrix(size, std::vector<T>(size));
#else
        std::vector<T> matrix(size*size);
#endif

#ifdef _ASMATRIX
        for(size_t index = 0; index < size; ++index) {
            matrix[index][index] = 1;
        }
#else
        for(size_t index = 0; index < size; ++index) {
            matrix[index + index*size + index] = 1;
        }
#endif

#ifdef _ASMATRIX
        return matrix;
#else
        return Matrix(size, size, matrix);
#endif
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

#ifdef _ASMATRIX
        std::vector<std::vector<T>> matrix(range.rows, std::vector<T>(range.columns));
#else
        std::vector<T> matrix(range.rows*range.columns);
#endif

        std::mt19937 mt(std::random_device{}());
        std::uniform_real_distribution<> real_dist(range.from, range.to);
        const auto gen = std::bind(std::ref(real_dist), std::ref(mt));

#ifdef _ASMATRIX
        for(auto & row : matrix) {
            std::generate(row.begin(), row.end(), [gen]() -> auto {
                return gen();
            });
        }
#else
        std::generate(matrix.begin(), matrix.end(), [gen]() -> auto {
            return gen();
        });
#endif

#ifdef _ASMATRIX
        return matrix;
#else
        return Matrix(range.rows, range.columns, matrix);
#endif
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
#ifdef _ASMATRIX
            {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
#else
            4, 4, {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}
#endif
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
