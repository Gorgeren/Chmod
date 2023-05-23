#include <iostream>
#include <vector>
#include <iomanip>
#include "lib/ChmodLib.h"
int main() {
    using namespace std;
    // Решаем задачу вида Ax = b
    ChmodLib::matrix A({{-1, -3, -4,  0},
                        { 3,  7, -8,  3},
                        { 1, -6,  2,  5},
                        {-8, -4, -1, -1}});

    vector<double> b = {-1, -3, -4, 0};
    auto x = ChmodLib::LU::solve_SLAU(A, b);
    ChmodLib::matrix answer(x);
    for(int i = 0; i < answer.row(); i++) {
        cout << "x" << i+1 << " = " << answer[0][i] << '\n';
    }
    answer = ChmodLib::transpose(answer);
    cout << "check answer: A*x\n";
    cout << A * answer;
    ChmodLib::matrix inversedA = ChmodLib::inverse(A);
    cout << "inversed A\n";
    cout << inversedA;
    cout << "Check that the inverse matrix is correct: A * A^-1\n";
    cout << A * inversedA;
    cout << "determinant of matrix A\n" <<ChmodLib::determinant(A) << endl;
    
}
