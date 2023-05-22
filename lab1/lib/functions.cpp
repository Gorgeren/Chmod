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
std::vector<double> SWEEP::solve_SLAU(const std::vector<double> &a, const std::vector<double>& b, const std::vector<double> &c, const std::vector<double>& d) {
    if(a.size() != b.size() && a.size() != c.size() && a.size() != d.size() && a.size() < 3) throw std::invalid_argument("it is not tridiagonal matrix");
    int n = a.size();
    std::vector<double> P(n);
    std::vector<double> Q(n);
    std::vector<double> ans(n);
    P[0] = -c[0] / b[0];
    Q[0] =  d[0] / b[0];
    for(int i = 1; i < n; i++) {
        P[i] = -c[i]/(b[i] + a[i] * P[i - 1]);
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
    }
    ans[n - 1] = Q[n - 1];

    for(int i = n - 2; i >= 0; i--) {
        ans[i] = P[i] * ans[i + 1] + Q[i];
    }
    return ans;
}
matrix SWEEP::makeMatrix(const std::vector<double> &a, const std::vector<double>& b, const std::vector<double> &c) {
    if(a.size() != b.size() && a.size() != c.size() && a.size() < 3) throw std::invalid_argument("it is not tridiagonal matrix");
    matrix res(a.size(), a.size());
    for(int i = 0; i < (int)a.size(); i++) {
        if(i != 0 && i != (int)a.size() - 1) {
            res[i][i - 1] = a[i];
            res[i][i] = b[i];
            res[i][i + 1] = c[i];
        }
        else if(i == 0) {
            res[i][i] = b[i];
            res[i][i + 1] = c[i];
        }
        else if(i == (int)a.size() - 1 ) {
            res[i][i] = b[i];
            res[i][i - 1] = a[i];
        }
    }
    return res;
}
std::pair<matrix, int> ITER::solve_SLAU_simple(ChmodLib::matrix &A, ChmodLib::matrix &B, double eps) {
    int iterations = 0;
    int n = A.col();
    matrix beta(n, 1);
    matrix alpha(n, n);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i != j) {
                if(A[i][i] == 0) {
                    std::cout << "did not have time to make a method with permutations. This is division by zero.\n";
                    exit(0);
                }
                alpha[i][j] = - A[i][j] / A[i][i];
            } else alpha[i][j] = 0;
        }
    }

    for(int i = 0; i < n; i++) {
        beta[i][0] = B[i][0] / A[i][i];
    }
    double g = 0;
    for(int i = 0; i < n; i++) {
        double h = 0;
        for(int j = 0; j < n; j++) {
            h += std::abs(alpha[i][j]);
        }
        if(h > g) g = h;
    }
    if(g >= 1) {
        std::cout << "The convergence condition is violated\n";
        exit(0);
    }
    double q = alpha.normC();
    double normB = beta.normC();
    double f = ITER::eps(q, normB);
    matrix xprev = beta;
    matrix x;
    // std::cout << "alpha\n";
    // std::cout << alpha;
    // std::cout << "beta\n";
    // std::cout << beta << '\n';
    while(f > eps) {
        iterations++;
        x = beta + alpha * xprev;
        // std::cout << x << '\n';
        matrix delta = x - xprev;
        normB = delta.normC();
        f = ITER::eps(q, normB);
        xprev = x;
    }
    return {x, iterations};
}
std::pair<matrix, int> ITER::solve_SLAU_seidel(matrix &A, matrix &B, double eps) {
    int iterations = 0;
    int n = A.col();
    matrix beta(n, 1);
    matrix alpha(n, n);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i != j) {
                if(A[i][i] == 0) {
                    std::cout << "did not have time to make a method with permutations. This is division by zero.\n";
                    exit(0);
                }
                alpha[i][j] = - A[i][j] / A[i][i];
            } else alpha[i][j] = 0;
        }
    }

    for(int i = 0; i < n; i++) {
        beta[i][0] = B[i][0] / A[i][i];
    }
    double g = 0;
    for(int i = 0; i < n; i++) {
        double h = 0;
        for(int j = 0; j < n; j++) {
            h += std::abs(alpha[i][j]);
        }
        if(h > g) g = h;
    }
    if(g >= 1) {
        std::cout << "The convergence condition is violated\n";
        exit(0);
    }
    double q = alpha.normC();
    double normB = beta.normC();
    double f = ITER::eps(q, normB);
    matrix xprev = beta;
    matrix x = beta;
    while(f > eps) {
        iterations++;
        for(int i = 0; i < n; i++) {
            x[i][0] = beta[i][0];
            for(int j = 0; j < n; j++) {
                x[i][0] += alpha[i][j] * x[j][0];
            }
        }
        matrix delta = x - xprev;
        normB = delta.normC();
        f = ITER::eps(q, normB);
        xprev = x;
    }
    return {x, iterations};
}
double ITER::eps(double q, double normB) {
    return (std::abs(q)/(1 - std::abs(q))) * normB;
}
}