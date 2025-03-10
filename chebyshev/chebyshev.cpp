#include <iostream>
#include <vector>
#include <cmath>
#include <locale> // For UTF-8 support

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to solve a system of linear equations using Chebyshev acceleration
void chebyshev_acceleration(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.size();
    std::vector<double> x_prev = x;

    for (int k = 1; k <= max_iter; ++k) {
        // Compute the next approximation
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

        // Chebyshev acceleration
        double omega = 1.0 / (1.0 - 0.25 * std::pow(M_PI / (2 * max_iter + 1), 2));
        for (int i = 0; i < n; ++i) {
            x[i] = x_prev[i] + omega * (x_next[i] - x_prev[i]);
        }

        // Check for convergence
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += std::abs(x[i] - x_prev[i]);
        }
        if (error < tol) {
            std::cout << "Chebyshev acceleration converged in " << k << " iterations." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::cout << "Chebyshev acceleration did not converge in " << max_iter << " iterations." << std::endl;
}

// Function to solve a system of linear equations using the Jacobi method
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
            std::cout << "Jacobi method converged in " << k + 1 << " iterations." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::cout << "Jacobi method did not converge in " << max_iter << " iterations." << std::endl;
}

// Function to solve a system of linear equations using the Gauss-Seidel method
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
            std::cout << "Gauss-Seidel method converged in " << k + 1 << " iterations." << std::endl;
            return;
        }

        x_prev = x;
    }

    std::cout << "Gauss-Seidel method did not converge in " << max_iter << " iterations." << std::endl;
}

int main() {
    // Set locale for UTF-8 support (if needed)
    std::locale::global(std::locale("en_US.UTF-8"));
    std::cout.imbue(std::locale());

    // Example system of linear equations
    std::vector<std::vector<double>> A = { {4, 1, 1}, {1, 4, 1}, {1, 1, 4} };
    std::vector<double> b = { 6, 6, 6 };
    std::vector<double> x = { 0, 0, 0 };

    int max_iter = 1000;
    double tol = 1e-6;

    // Solve using Chebyshev acceleration
    chebyshev_acceleration(A, b, x, max_iter, tol);

    // Solve using Jacobi method
    x = { 0, 0, 0 };
    jacobi_method(A, b, x, max_iter, tol);

    // Solve using Gauss-Seidel method
    x = { 0, 0, 0 };
    gauss_seidel_method(A, b, x, max_iter, tol);

    return 0;
}
