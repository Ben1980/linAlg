#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>

#include "matrix.h"

int main(int argc, char **argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    int res = context.run();

    Matrix<double, 3, 3> a = {{}};

    if (context.shouldExit()) {
        return res;
    }
}

