#ifndef LINALG_UTILS_H
#define LINALG_UTILS_H

//#include <fmt/format.h>
//#include <vector>
//#include <algorithm>
//#include "matrix.h"

namespace Utils {
    template<typename T>
    void PrintMatrix(const Matrix<T> &m) {
        for(const auto &row : m()) {
            //fmt::print("| {:^5} |\n", fmt::join(row, ""));
        }
    }
}

#endif //LINALG_UTILS_H
