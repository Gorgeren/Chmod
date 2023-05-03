#include <iostream>
#include <vector>
#include <iomanip>
#include "lib/ChmodLib.h"
int main() {
    using namespace std;
    ChmodLib::matrix A({{1000, 2 ,3},
                        {2, 100, 5},
                        {400, 2, 3}
                       });
    cout << "matrix A\n";
    cout << A;
    auto [L, U, P, chet] = ChmodLib::makeLU(A);
    cout << "matrix L\n";
    cout << L;
    cout << "matrix U\n";
    cout << U;
    cout << "determinant of matrix A == ";
    cout << ChmodLib::LU::determinant(U, chet) << endl;
    cout << ChmodLib::permutate(L*U, P);
    
}
