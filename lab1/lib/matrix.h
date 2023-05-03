#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
namespace ChmodLib {
    class matrix {
    public:
        matrix() = default;
        matrix(int n, int m);
        matrix(std::vector<std::vector<double>> &t);
        matrix(const std::initializer_list<std::initializer_list<double>> &list);

        matrix(std::pair<int, int> size): n(size.first), m(size.second) {};
        int col() const {return n;}
        int row() const {return m;}
        void print();
        std::pair<int, int> size() const;
        matrix operator*(const matrix& p) const;
        matrix operator-(const int num) const;
        matrix operator*(const int num) const;
        std::vector<double>& operator[](const int index);
        const std::vector<double>& operator[](int index) const;
    private:
        std::vector<std::vector<double>> mat;
        int n, m;
        friend matrix operator*(const int num, const matrix& a);
        friend std::ostream& operator<<(std::ostream& os, const matrix& mat);
    };
}