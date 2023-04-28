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
std::pair<matrix, matrix> LU(const matrix& mat) {
    matrix L(mat.size().first, mat.size().second);
    for(int i = 0; i < std::min(L.col(), L.row()); i++) {
        L[i][i] = 1;
    }
    matrix tmp = mat;
    matrix check;
    matrix P = make_identity(std::min(L.col(), L.row()));
    for(int j = 0; j < L.row(); j++) {
        int index = find_index(tmp, j, j);
        if(index != -1) {
            std::cout << "first check" << std::endl << L *tmp;
            std::cout << "matrix P" << std::endl << P;
            std::swap(tmp[index], tmp[j]);
            // std::swap(L[index], L[j]);
            std::swap(P[index], P[j]);
            std::cout << "after swaping"  << std::endl<< P;
        } else throw std::invalid_argument("matrix is degenerate");
        check = make_identity(L.col());
        for(int i = j; i < L.col() - 1; i++) {
            double m = tmp[i + 1][j] / tmp[j][j];
            double tmp_i_1_j = tmp[i + 1][j];
            double tmp_j_j = tmp[j][j];
            std::cout << "befor" << std::endl;
            std::cout <<tmp;
            L[i + 1][j] = m;
            check[i + 1][j] = m;
            for(int k = j; k < L.row(); k++) {
                tmp[i + 1][k] -= tmp[j][k] * m;
            }
            std::cout << "after" << std::endl;
            std::cout << tmp;
            std::cout << "matrix L" << std::endl;
            std::cout << L;
        }
    }
    return {P * L,  tmp};

}
std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b) {

}

}