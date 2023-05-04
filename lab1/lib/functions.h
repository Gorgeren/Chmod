#pragma once
#include "matrix.h"
#include <algorithm>
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
        std::vector<double> solve_SLAU(std::vector<double> &a, std::vector<double>& b, std::vector<double> &c, std::vector<double>& d);
    }

}