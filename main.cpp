#define DOCTEST_CONFIG_IMPLEMENT

#include <doctest/doctest.h>
#include <fmt/format.h>
#include <chrono>
#include <array>
#include "matrix.h"
#include "matrixFactory.h"
#include "luDecomposition.h"
#include "pivotLUDecomposition.h"

int main(int argc, char **argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run();

    Matrix<double> a = {
        3, 3, (std::array<double, 9>{3, 2, 1, 1, 0, 2, 1, 1, 1}).data()
    };
    Matrix<int> identity = MatrixFactory::IdentityMatrix<int>(4);

    fmt::print("\n");

    const size_t maxNumberOfMatrixElements = 1024;
#ifdef _ASMATRIX
    fmt::print("Performance test of matrix impelemntation with two dimensional vector\n");
    fmt::print("---------------------------------------------------------------------\n");
    for(size_t i = 2; i <= maxNumberOfMatrixElements; i *= 2) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const double elapsed= std::chrono::duration<double>(end - start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {:.10f}s\n", i,i, elapsed);
    }
#elif _ASARRAY
    fmt::print("Performance test of matrix impelemntation with raw array\n");
    fmt::print("--------------------------------------------------------\n");
    for(size_t i = 2; i <= maxNumberOfMatrixElements; i *= 2) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const double elapsed= std::chrono::duration<double>(end - start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {:.10f}s\n", i,i, elapsed);
    }
#else
    fmt::print("Performance test of matrix impelemntation with one dimensional vector\n");
    fmt::print("---------------------------------------------------------------------\n");
    for(size_t i = 2; i <= maxNumberOfMatrixElements; i *= 2) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const double elapsed= std::chrono::duration<double>(end - start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {:.10f}s\n", i,i, elapsed);
    }
#endif

    fmt::print("\n");

    LUDecomposition::Decomposition LU = LUDecomposition::Decompose(a);

    PivotLUDecomposition::Decomposition pivotLU = PivotLUDecomposition::Decompose(a);

    if (context.shouldExit()) {
        return res;
    }
}
