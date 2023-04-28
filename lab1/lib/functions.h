#pragma once
#include "matrix.h"
namespace ChmodLib {
    matrix transpose(const matrix& mat);
    std::pair<matrix, matrix> LU(const matrix& mat);
    std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b);
}