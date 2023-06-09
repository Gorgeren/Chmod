#pragma once
#include "matrix.h"
#include <algorithm>
#include <cmath>
namespace ChmodLib {
    // транспонирование матриц
    matrix transpose(const matrix& mat);
    // возвращает L, U, P.
    // L - нижнетреугольная матрица
    // U - верхнетреугольная
    // P - вектор, в котором хранится информация о перестановках, сделанных при разложении
    // chet - информация о четности количества перестановок строк
    std::tuple<matrix, matrix, std::vector<int>, bool> makeLU(const matrix& mat); 
    // Решает СЛАУ с использованием LU разложения
    namespace LU {
        std::vector<double> solve_SLAU(const matrix& A, const std::vector<double>& b);
    }
    // Делает обратные перестановки, которые сделаны на векторе P
    template <typename T>
    T permutate(T arr, std::vector<int> P) {
        for(int i = 0; i < P.size(); i++) {
            while(P[i] != i) {
                std::swap(arr[i], arr[P[i]]);
                std::swap(P[i], P[P[i]]);
            }
        }
        return arr;
    }
    template <typename T>
    //Перемешивает элементы, согласно индексам из P
    T permutate2(const T& arr, std::vector<int> P) {
        T smth = arr;
        for(int i = 0; i < (int)P.size(); i++) {
            smth[i] = arr[P[i]];
        }
        return smth;
    }
    namespace LU{
        // Возвращает детерминант матрицы, используя матрицу U
        // и информацию о четности количества перестановок
        double determinant(const matrix& U, bool chet);
    }
    // Возвращает детерминант матрицы
    double determinant(const matrix& A);
    // Обратная матрица
    matrix inverse(const matrix& mat);
    namespace SWEEP { // метод прогонки
        //решение СЛАУ трехдиогональной матрицы
        std::vector<double> solve_SLAU(const std::vector<double> &a, const std::vector<double>& b, const std::vector<double> &c, const std::vector<double>& d);
        //перевод трехдиогональной матрицы, задаваемой тремя векторами в обычную матрицу
        matrix makeMatrix(const std::vector<double> &a, const std::vector<double>& b, const std::vector<double> &c);
    }
    namespace ITER { // Итерационные методы
        //Нахождение eps(k)
        double eps(double q, double beta);
        //Решение СЛАУ методом итераций Зейделя
        std::pair<matrix, int> solve_SLAU_seidel(matrix &A, matrix &B, double eps);
        //Решение СЛАУ методом простых итераций
        std::pair<matrix, int> solve_SLAU_simple(matrix &A, matrix &B, double eps);
        //Возвращает матрицу U для метода Якоби
        // matrix U_matrix(int n, int i, int j, double cos, double sin);
        //Метод вращения Якоби. Возвращает вектор сз, матрицу св и кол-во итераций
        std::tuple<matrix, matrix, int> YakobiMethod(const matrix& A, double eps); 

    }
}