#pragma once
#include "matrix.h"
namespace ChmodLib {
    matrix transpose(const matrix& mat);
    
    std::tuple<matrix, matrix, std::vector<std::pair<int, int>>> LU(const matrix& mat); //returns L, U, P
    std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b);
}