#pragma once
#include "matrix.h"
#include <algorithm>
namespace ChmodLib {
    matrix transpose(const matrix& mat);
    std::tuple<matrix, matrix, std::vector<int>, bool> makeLU(const matrix& mat); //returns L, U, P
    std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b);
    
    template <typename T>
    T permutate(T arr, std::vector<int> P) {
        for(int i = 0; i < P.size(); i++) {
            while(P[i] != i) {
                std::swap(arr[i], arr[P[i]]);
                std::swap(P[i], P[P[i]]);
            }
        }
        return arr;
    }
    namespace LU{
        double determinant(const matrix& U, bool chet);
    }
}