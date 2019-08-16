#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>

#include "matrix.h"

int main(int argc, char **argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run();

    Matrix<double, 2, 3> a = {{{3,2,1},{1,0,2}}};

    fmt::print("\n");
    Helper::PrintMatrix(a);

    if (context.shouldExit()) {
        return res;
    }
}

