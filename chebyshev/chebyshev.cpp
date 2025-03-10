#include <iostream>
#include <vector>
#include <cmath>
#include <locale> // Для поддержки UTF-8

// Определение M_PI, если оно не определено
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Функция для решения системы линейных уравнений методом простых итераций с чебышёвским ускорением
void chebyshev_acceleration(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.size();
    std::vector<double> x_prev = x;

    for (int k = 1; k <= max_iter; ++k) {
        // Вычисление следующего приближения
        std::vector<double> x_next(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x_prev[j];
                }
            }
            x_next[i] = (b[i] - sum) / A[i][i];
        }

        // Чебышёвское ускорение
        double omega = 1.0 / (1.0 - 0.25 * std::pow(M_PI / (2 * max_iter + 1), 2));
        for (int i = 0; i < n; ++i) {
            x[i] = x_prev[i] + omega * (x_next[i] - x_prev[i]);
        }

        // Проверка на сходимость
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x[i] - x_prev[i]);
        }
        if (error < tol) {
            std::wcout << L"Чебышёвское ускорение сошлось за " << k << L" итераций." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"Чебышёвское ускорение не сошлось за " << max_iter << L" итераций." << std::endl;
}

// Функция для решения системы линейных уравнений методом Якоби
void jacobi_method(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.size();
    std::vector<double> x_prev = x;

    for (int k = 0; k < max_iter; ++k) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x_prev[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x[i] - x_prev[i]);
        }
        if (error < tol) {
            std::wcout << L"Метод Якоби сошёлся за " << k + 1 << L" итераций." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"Метод Якоби не сошёлся за " << max_iter << L" итераций." << std::endl;
}

// Функция для решения системы линейных уравнений методом Гаусса-Зейделя
void gauss_seidel_method(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.size();
    std::vector<double> x_prev = x;

    for (int k = 0; k < max_iter; ++k) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x[i] - x_prev[i]);
        }
        if (error < tol) {
            std::wcout << L"Метод Гаусса-Зейделя сошёлся за " << k + 1 << L" итераций." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"Метод Гаусса-Зейделя не сошёлся за " << max_iter << L" итераций." << std::endl;
}

int main() {
    // Установка локали для поддержки UTF-8
    std::locale::global(std::locale("en_US.UTF-8"));
    std::wcout.imbue(std::locale());

    // Пример системы линейных уравнений
    std::vector<std::vector<double>> A = { {4, 1, 1}, {1, 4, 1}, {1, 1, 4} };
    std::vector<double> b = { 6, 6, 6 };
    std::vector<double> x = { 0, 0, 0 };

    int max_iter = 1000;
    double tol = 1e-6;

    // Решение с использованием чебышёвского ускорения
    chebyshev_acceleration(A, b, x, max_iter, tol);

    // Решение с использованием метода Якоби
    x = { 0, 0, 0 };
    jacobi_method(A, b, x, max_iter, tol);

    // Решение с использованием метода Гаусса-Зейделя
    x = { 0, 0, 0 };
    gauss_seidel_method(A, b, x, max_iter, tol);

    return 0;
}