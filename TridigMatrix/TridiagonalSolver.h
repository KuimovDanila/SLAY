#ifndef TRIDIAGONAL_SOLVER_H
#define TRIDIAGONAL_SOLVER_H

#include <vector>

class TridiagonalSolver {
public:
    TridiagonalSolver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d);
    std::vector<double> solve();

private:
    std::vector<double> a; // Нижняя диагональ (под главной)
    std::vector<double> b; // Главная диагональ
    std::vector<double> c; // Верхняя диагональ (над главной)
    std::vector<double> d; // Правая часть уравнения
};

#endif // TRIDIAGONAL_SOLVER_H