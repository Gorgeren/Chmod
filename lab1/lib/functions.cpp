#include "functions.h"
#include <algorithm>
namespace ChmodLib {
matrix transpose(const matrix& matr) {
    std::pair<int, int> size = matr.size();
    matrix res(size.second, size.first);

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
    for(int i = 0; i < (int)P.size(); i++) P[i] = i;
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
double determinant(const matrix& A) {
    auto [L, U, P, chet] = makeLU(A);
    if(U.col() != U.row()) throw std::invalid_argument("The determinant of a non-square matrix does not exist");
    double res = 1;
    for(int i = 0; i < U.row(); i++) {
        res *= U[i][i];
    }
    return chet ? -1 * res: res;
}
std::vector<double> LU::solve_SLAU(const matrix& A, const std::vector<double>& b) {
    auto [L, U, P, chet] = makeLU(A);
    std::vector<double> b_ = permutate2(b, P);

    // Решаем задачу Lz = b, z = Ux.
    std::vector<double> z(L.row());
    z[0] = b_[0];
    for(int i = 1; i < L.col(); i++) {
        z[i] = b_[i];
        for(int j = 0; j < i; j++) {
            z[i] -= z[j] * L[i][j];
        }
    }

    // Решаем задачу Ux = z
    std::vector<double> x(L.col());
    int last = (int)x.size() - 1;
    x[last] = z[last] / U[last][last];
    for(int i = last - 1; i >= 0; i--) {
        x[i] = z[i];
        for(int j = last; j > i; j--) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}

matrix inverse(const matrix& mat) {
    matrix tmp = make_identity(mat.col());
    std::vector<double> tmp_b(mat.col());
    matrix ans(mat.col(), mat.col());
    for(int i = 0; i < (int)tmp_b.size(); i++) {
        if(i != 0) {
            tmp_b[i - 1] = 0;
            tmp_b[i] = 1;
        } tmp_b[i] = 1;
        std::vector<double> tmp = LU::solve_SLAU(mat, tmp_b);
        ans[i] = tmp;
    }
    ans = transpose(ans);
    return ans;
}
std::vector<double> SWEEP::solve_SLAU(std::vector<double> &a, std::vector<double>& b, std::vector<double> &c, std::vector<double>& d) {
    if(a.size() != b.size() && a.size() != c.size() && a.size() != d.size() && a.size() < 3) throw std::invalid_argument("it is not tridiagonal matrix");
    std::vector<double> P(a.size());
    std::vector<double> Q(a.size());
    std::vector<double> ans(a.size());
    P[0] = -c[0] / b[0];
    Q[0] =  d[0] / b[0];
    for(int i = 1; i < (int)a.size() - 1; i++) {
        P[i] = -c[i]/(b[i] + a[i] * P[i - 1]);
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
    }
    ans[a.size() - 1] = Q[a.size() - 1];
    for(int i = a.size() - 2; i > 0; i--) {
        ans[i] = P[i] * ans[i + 1] + Q[i];
    }
    return ans;
}
        
}