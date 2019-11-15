#ifndef LINALG_UTILS_H
#define LINALG_UTILS_H

#include <fmt/format.h>
#include "matrix.h"

namespace Utils {
    template<typename T>
    void PrintMatrix(const Matrix<T> &m) {
        for(size_t row = 0; row < m.rows(); ++row) {
            for(size_t column = 0; column < m.columns(); ++column) {
                fmt::print("{} ", m(row, column));
            }
            fmt::print("\n");
        }
        fmt::print("\n");
    }
}

#endif //LINALG_UTILS_H
