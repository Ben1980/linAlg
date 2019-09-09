#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>

#include <chrono>
#include "matrix.h"
#include "utils.h"
#include "matrixFactory.h"
#include "luDecomposition.h"
#include "gauss.h"

int main(int argc, char **argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run();

    Matrix<double> a = {{{3,2,1},{1,0,2}}};
    Matrix<int> identity = MatrixFactory::IdentityMatrix<int>(4);

    fmt::print("\n");

    for(size_t i = 2; i <= 256; i *= i) {
        Matrix<double> A = MatrixFactory::RandomMatrix<double>({i, i, -100., 100.});

        const auto start = std::chrono::system_clock::now();
        Matrix C = A * A;
        const auto end = std::chrono::system_clock::now();
        const int elapsed= std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
        fmt::print("Elapsed time for C = A * A, {}x{}: {}ns\n", i,i, elapsed);
    }
    fmt::print("\n");

    fmt::print("LU-Decomposition with pivoting");
    LUDecomposition::Decomposition test = LUDecomposition::Decompose(a);

    Matrix<double> b = Gauss::Solve<double>(a, {-10,-19,-11});

    if (context.shouldExit()) {
        return res;
    }
}
