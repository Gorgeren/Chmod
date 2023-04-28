#include <iostream>
#include <vector>
#include <iomanip>
#include "lib/ChmodLib.h"
int main() {
    using namespace std;
    ChmodLib::matrix test({{1, 2 ,3},
                            {2, 0, 5},
                            {4, 2, 3}
                          });
    cout << test;

    auto [L, U] = ChmodLib::LU(test);
    cout << "cheking\n";
    cout << L * U;
    cout << "matrix L\n";
    std::cout << L;
    cout << "matrix U\n";
    std::cout << U;
    // ChmodLib::matrix first({{1, 0, 0, 0},
    //                         {2, 1, 0, 0},
    //                         {3, 0, 1, 0},
    //                         {0, 0, 0, 1}});
    // ChmodLib::matrix second({{1, 0, 0, 0},
    //                          {0, 1, 0, 0},
    //                          {0, 2, 1, 0},
    //                          {0, 2, 0, 1}});

}
