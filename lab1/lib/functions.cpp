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
    bool flag = 1;
    for(int i = _i; i < mat.size().first; i++) {
        flag = 0;
        if(std::abs(mat[i][_j]) > std::abs(max)) {
            max = std::abs(mat[i][_j]);
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
std::tuple<matrix, matrix, std::vector<int>, bool> makeLU(const matrix& mat) {
    bool chet = 0;
    std::vector<int> P(mat.size().first);
    for(int i = 0; i < P.size(); i++) P[i] = i;
    matrix L(mat.size().first, mat.size().second);
    matrix U = mat;
    for(int i = 0; i < std::min(L.col(), L.row()); i++) {
        L[i][i] = 1;
    }
    for(int j = 0; j < L.row(); j++) {
        int index = find_index(U, j, j);
        if(index != -1) {
            if(index != j) {
                chet = !chet;
                std::swap(U[index], U[j]);
                std::swap(P[index], P[j]);
                std::swap(L[index], L[j]);
                swap_column(L, index, j);
            }
        } else throw std::invalid_argument("matrix is degenerate");
        for(int i = j; i < L.col() - 1; i++) {
            double m = U[i + 1][j] / U[j][j];
            L[i + 1][j] = m;
            for(int k = j; k < L.row(); k++) {
                U[i + 1][k] -= U[j][k] * m;
            }
        }
        std::cout << "matrix U\n";
        std::cout << U;
        std::cout << "matrix L\n";
        std::cout << L;
    }
    return {L, U, P, chet};
}
double LU::determinant(const matrix& U, bool chet) {
    if(U.col() != U.row()) throw std::invalid_argument("The determinant of a non-square matrix does not exist");
    double res = 1;
    for(int i = 0; i < U.row(); i++) {
        res *= U[i][i];
    }
    return chet ? -1 * res: res;
}
std::vector<double> solve_SLAU(const matrix &L, const matrix &U, const std::vector<double>& b) {

}


}