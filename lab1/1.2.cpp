#include <iostream>
#include <vector>
#include "lib/ChmodLib.h"
using namespace std;
int main() {
    vector<double> a = {0, 7, -9, 7, -4};
    vector<double> b = {-1, -17, 19, -20, 12};
    vector<double> c = {-1, -8, 8, 4, 0};
    vector<double> d = {-4, 132, -59, -193, -40};
    vector<double> x = ChmodLib::SWEEP::solve_SLAU(a, b, c, d);
    ChmodLib::matrix A = {{-1, -1,  0, 0, 0},
                          {7, -17, -8, 0, 0},
                          {0,  -9, 19, 8, 9},
                          {0,  0, 7, -20, 4},
                          {0,  0, 0, -4, 12}};
    cout << "test\n";
    cout << A;
    cout << "first ans\n";
    for(int i = 0; i < (int)x.size(); i++) {
        cout << x[i] << ' ';
    }
    cout << '\n';
    ChmodLib::matrix ch(x);
    ch = ChmodLib::transpose(ch);
    cout << "check ans\n";
    cout << A * ch;
    x = ChmodLib::LU::solve_SLAU(A, d);
    cout << "lu ans\n";
    for(int i: x) cout << i << ' ';
    cout << '\n';
    ChmodLib::matrix l(x);
    l = ChmodLib::transpose(l);
    cout <<"lu check\n";
    cout << A * l;

    ChmodLib::matrix ttt = ChmodLib::SWEEP::makeMatrix(a,b,c);
    cout << "matrix A\n";
    cout << A;
    cout << "matrix from vectors\n";
    cout << ttt;
    cout << "_________\n";
    cout << ttt * l;
    // b1 c1 0  0  0  0  d1 -> x1 = (d1 - c1x2) / b1 
    // a2 b2 c2 0  0  0  d2 -> a2((d1 - c1) / b1) + b2x + c2
    // 0  a3 b3 c3 0  0  d3
    // 0  0  a4 b4 c4 0  d4
    // 0  0  0  a5 b5 c5 d5
    // 0  0  0  0  а6 b6 d6


    
    // x1 = -c1/b1 + d1/b1 -> P1 = -c1/b1, Q1 = d1/b1


}   