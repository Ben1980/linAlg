#define DOCTEST_CONFIG_IMPLEMENT

#include <doctest/doctest.h>
#include <fmt/format.h>
#include <chrono>
#include <array>
#include "matrix.h"
#include "matrixFactory.h"
#include "luDecomposition.h"
#include "gauss.h"

int main(int argc, char **argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run();

    Matrix<double> a = {
#ifdef _ASMATRIX
        {{3,2,1},{1,0,2}}
#elif _ASARRAY
        2, 3, (std::array<double, 6>{3, 2, 1, 1, 0, 2}).data()
#else
        2, 3, {3, 2, 1, 1, 0, 2}
#endif
    };
    Matrix<int> identity = MatrixFactory::IdentityMatrix<int>(4);

    fmt::print("\n");

#ifdef _ASMATRIX
    fmt::print("Performance test of matrix impelemntation with two dimensional vector\n");
    fmt::print("---------------------------------------------------------------------\n");
    for(size_t i = 2; i <= 256; i *= i) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const int elapsed= std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {}ns\n", i,i, elapsed);
    }
#elif _ASARRAY
    fmt::print("Performance test of matrix impelemntation with raw array\n");
    fmt::print("--------------------------------------------------------\n");
    for(size_t i = 2; i <= 256; i *= i) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const int elapsed= std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {}ns\n", i,i, elapsed);
    }
#else
    fmt::print("Performance test of matrix impelemntation with one dimensional vector\n");
    fmt::print("---------------------------------------------------------------------\n");
    for(size_t i = 2; i <= 256; i *= i) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const int elapsed= std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {}ns\n", i,i, elapsed);
    }
#endif

    fmt::print("\n");

    fmt::print("LU-Decomposition with pivoting");
    LUDecomposition::Decomposition test = LUDecomposition::Decompose(a);

    Gauss::Decomposition<double> b = Gauss::Decompose(a);

    if (context.shouldExit()) {
        return res;
    }
}
