#include <iostream>
#include <vector>
#include "lib/ChmodLib.h"
using namespace std;

int main() {
    // ChmodLib::matrix A{{-7, -5, -9},
    //                    {-5,  5,  2},
    //                    { 9,  2,  9}};
    ChmodLib::matrix A{{4, 2, 1},
                       {2, 5, 3},
                       {1, 3, 6}};
    double eps;
    cout << "Enter eps\n";
    cin >> eps;
    double tmp = -6;
    cout << std::atan(tmp) << '\n';
    auto [sobstv, vectors, iter] = ChmodLib::ITER::YakobiMethod(A, eps);
    cout << sobstv << '\n';
    cout << vectors << '\n';
    cout << "iterations == " << iter << '\n';
    

}