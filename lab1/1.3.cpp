#include <iostream>
#include <vector>
#include "lib/ChmodLib.h"
using namespace std;
int main() {
    ChmodLib::matrix A({{-22, -2, -6, 6},
                        {3,  -17, -3, 7},
                        {2,    6, -17, 5},
                        {-1,  -8,  8, 23}});

    vector<double> b = {96, -26, 35, -234};
    ChmodLib::matrix B(b);
    B = ChmodLib::transpose(B);
    auto x = ChmodLib::LU::solve_SLAU(A, b);
    ChmodLib::matrix X(x);
    cout << X;
    double eps;
    cout << "ITER method\n";
    cout << "Enter eps:";
    cin >> eps;
    auto [ans, iter] = ChmodLib::ITER::solve_SLAU_simple(A, B, eps);
    cout << ChmodLib::transpose(ans);
    std::cout << "iterations: " << iter << '\n';
    cout << "Seidel method\n";
    auto [ans2, iter2] = ChmodLib::ITER::solve_SLAU_seidel(A, B, eps);
    cout << ChmodLib::transpose(ans2);
    std::cout << "iterations: " << iter2 << '\n';

}