#include "functions.h"
#include <algorithm>
namespace ChmodLib {
matrix transpose(const matrix& matr) {
    std::pair<int, int> size = matr.size();
    matrix res(size.first, size.second);

    for(int i = 0; i < size.second; i++) {
        for(int j = 0; j < size.first; j++) {
            res[i][j] = matr[j][i];
        }
    }
    return res;  
}
int find_index(const matrix& mat, int _i, int _j) {
    double max = 0;
    int index = 0;
    for(int i = _i; i < mat.size().first; i++) {
        if(mat[i][_j] > std::abs(max)) {
            max = mat[i][_j];
            index = i;
        }
    }
    return max == 0 ? -1 : index;
}
matrix make_identity(int size) {
    matrix res(size, size);
    for(int i = 0; i < size; i++) {
        res[i][i] = 1;
    }
    return res;
}
void swap_column(matrix& mat, int i, int j) {
    int size = mat.col();
    double tmp;
    for(int k = 0; k < size; k++) {
        tmp = mat[k][i];
        mat[k][i] = mat[k][j];
        mat[k][j] = tmp;
    }
}
std::tuple<matrix, matrix, std::vector<std::pair<int, int>>> LU(const matrix& mat) {
    std::vector<std::pair<int, int>> P;
    matrix L(mat.size().first, mat.size().second);
    matrix U = mat;
    for(int i = 0; i < std::min(L.col(), L.row()); i++) {
        L[i][i] = 1;
    }
    for(int j = 0; j < L.row(); j++) {
        int index = find_index(U, j, j);
        if(index != -1 && index != j) {
            std::swap(U[index], U[j]);
            std::swap(L[index], L[j]);
            swap_column(L, index, j);
            P.push_back({j, index});
        } else throw std::invalid_argument("matrix is degenerate");
        for(int i = j; i < L.col() - 1; i++) {
            double m = U[i + 1][j] / U[j][j];
            L[i + 1][j] = m;
            for(int k = j; k < L.row(); k++) {
                U[i + 1][k] -= U[j][k] * m;
            }
        }
    }
    return {L, U, P};
}
std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b) {

}

}