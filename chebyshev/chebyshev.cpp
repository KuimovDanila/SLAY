#include <iostream>
#include <vector>
#include <cmath>
#include <locale> // ��� ��������� UTF-8

// ����������� M_PI, ���� ��� �� ����������
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ������� ��� ������� ������� �������� ��������� ������� ������� �������� � ����������� ����������
void chebyshev_acceleration(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.size();
    std::vector<double> x_prev = x;

    for (int k = 1; k <= max_iter; ++k) {
        // ���������� ���������� �����������
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

        // ����������� ���������
        double omega = 1.0 / (1.0 - 0.25 * std::pow(M_PI / (2 * max_iter + 1), 2));
        for (int i = 0; i < n; ++i) {
            x[i] = x_prev[i] + omega * (x_next[i] - x_prev[i]);
        }

        // �������� �� ����������
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x[i] - x_prev[i]);
        }
        if (error < tol) {
            std::wcout << L"����������� ��������� ������� �� " << k << L" ��������." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"����������� ��������� �� ������� �� " << max_iter << L" ��������." << std::endl;
}

// ������� ��� ������� ������� �������� ��������� ������� �����
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
            std::wcout << L"����� ����� ������� �� " << k + 1 << L" ��������." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"����� ����� �� ������� �� " << max_iter << L" ��������." << std::endl;
}

// ������� ��� ������� ������� �������� ��������� ������� ������-�������
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
            std::wcout << L"����� ������-������� ������� �� " << k + 1 << L" ��������." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::wcout << L"����� ������-������� �� ������� �� " << max_iter << L" ��������." << std::endl;
}

int main() {
    // ��������� ������ ��� ��������� UTF-8
    std::locale::global(std::locale("en_US.UTF-8"));
    std::wcout.imbue(std::locale());

    // ������ ������� �������� ���������
    std::vector<std::vector<double>> A = { {4, 1, 1}, {1, 4, 1}, {1, 1, 4} };
    std::vector<double> b = { 6, 6, 6 };
    std::vector<double> x = { 0, 0, 0 };

    int max_iter = 1000;
    double tol = 1e-6;

    // ������� � �������������� ������������ ���������
    chebyshev_acceleration(A, b, x, max_iter, tol);

    // ������� � �������������� ������ �����
    x = { 0, 0, 0 };
    jacobi_method(A, b, x, max_iter, tol);

    // ������� � �������������� ������ ������-�������
    x = { 0, 0, 0 };
    gauss_seidel_method(A, b, x, max_iter, tol);

    return 0;
}