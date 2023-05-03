#include "matrix.h"
namespace ChmodLib {
void matrix::print() {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            std::cout << std::setw(4) << mat[i][j];
        }
        std::cout << '\n';
    }
}
std::pair<int, int> matrix::size() const {
    return {n, m};
}
matrix matrix::operator*(const matrix& p) const {
    if(m != p.n) throw std::invalid_argument("ERROR: matrixs size must be a x m and m x n");
    matrix res(n, p.m);
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < p.m; j++) {
            double sum = 0.0;
            for(int k = 0; k < m; k++) {
                sum += mat[i][k] * p.mat[k][j];
            }
            res.mat[i][j] = sum;

        }
    }
    return res;
}

matrix::matrix(const std::initializer_list<std::initializer_list<double>> &list): n(list.size()), m(0) {
    if(n > 0) {
        int size = list.begin()->size();
        m = size;
        for(const auto& row: list) {
            if(size != row.size()) throw std::invalid_argument("ERROR: it is not a matrix");
            mat.push_back(std::vector<double>(row));
        }
    }
}
matrix matrix::operator-(const int num) const {
    matrix res(n, m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res.mat[i][j] = mat[i][j] - num;
        }
    }
    return res;
}
std::vector<double>& matrix::operator[](const int index) {
    return mat[index];
}
const std::vector<double>& matrix::operator[](const int index) const{
    return mat[index];
}

matrix matrix::operator*(const int num) const {
    matrix res(n, m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res.mat[i][j] = mat[i][j] * num;
        }
    }
    return res;
}
matrix operator*(const int num, const matrix& a) {
    matrix res(a.n, a.m);
    for(int i = 0; i < a.n; i++) {
        for(int j = 0; j < a.m; j++) {
            res.mat[i][j] = a.mat[i][j] * num;
        }
    }
    return res;
}
std::ostream& operator<<(std::ostream& os, const matrix& mat) {
    for(int i = 0; i < mat.n; i++) {
        for(int j = 0; j < mat.m; j++) {
            os.precision(4);
            os << std::fixed <<std::setw(7) << mat.mat[i][j] << ' ';
        }
        os << '\n';
    }
    return os;
}
matrix::matrix(int n, int m) : n(n), m(m), mat(n, std::vector<double> (m)) {};
matrix::matrix(std::vector<std::vector<double>> &t) : n(t.size()), m(0),mat(t) {
    if(n > 0) {
        m = t[0].size();
        int tmp = m;
        for(int i = 0; i < t.size(); i++) {
            if(t[i].size() != tmp) throw std::invalid_argument("ERROR: it is not a matrix");
        }
    }
}
}
