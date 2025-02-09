#include <iostream>
#include "TridiagonalSolver.h"

int main() {
    std::vector<double> a = { 0, -1, -1, -1 }; // Нижняя диагональ
    std::vector<double> b = { 4, 4, 4, 4 };    // Главная диагональ
    std::vector<double> c = { -1, -1, -1, 0 }; // Верхняя диагональ
    std::vector<double> d = { 5, 5, 10, 23 };  // Правая часть

    TridiagonalSolver solver(a, b, c, d);
    std::vector<double> solution = solver.solve();

    std::cout << "Solution: ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}


