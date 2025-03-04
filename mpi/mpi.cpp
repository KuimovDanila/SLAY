#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

// Function to print a vector
void printVector(const vector<double>& vec) {
    for (const auto& val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

// Simple Iteration Method
vector<double> simpleIteration(const vector<vector<double>>& A, const vector<double>& b, double eps, int& iterations) {
    int n = A.size();
    vector<double> x(n, 0.0);
    vector<double> x_new(n, 0.0);
    double norm;

    iterations = 0;
    do {
        norm = 0.0;
        for (int i = 0; i < n; ++i) {
            x_new[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x_new[i] -= A[i][j] * x[j];
                }
            }
            x_new[i] /= A[i][i];
            norm += abs(x_new[i] - x[i]);
        }
        x = x_new;
        iterations++;
    } while (norm > eps);

    return x;
}

// Jacobi Method
vector<double> jacobiMethod(const vector<vector<double>>& A, const vector<double>& b, double eps, int& iterations) {
    int n = A.size();
    vector<double> x(n, 0.0);
    vector<double> x_new(n, 0.0);
    double norm;

    iterations = 0;
    do {
        norm = 0.0;
        for (int i = 0; i < n; ++i) {
            x_new[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x_new[i] -= A[i][j] * x[j];
                }
            }
            x_new[i] /= A[i][i];
            norm += abs(x_new[i] - x[i]);
        }
        x = x_new;
        iterations++;
    } while (norm > eps);

    return x;
}

// Gauss-Seidel Method
vector<double> gaussSeidelMethod(const vector<vector<double>>& A, const vector<double>& b, double eps, int& iterations) {
    int n = A.size();
    vector<double> x(n, 0.0);
    double norm;

    iterations = 0;
    do {
        norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double sum = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= A[i][j] * x[j];
                }
            }
            double x_new = sum / A[i][i];
            norm += abs(x_new - x[i]);
            x[i] = x_new;
        }
        iterations++;
    } while (norm > eps);

    return x;
}

int main() {
    // Example system of linear equations
    vector<vector<double>> A = { {10, -1, 2},
                                {-1, 11, -1},
                                {2, -1, 10} };
    vector<double> b = { 6, 25, -11 };
    double eps = 1e-6;

    int iterations;
    vector<double> x;

    // Simple Iteration Method
    auto start = high_resolution_clock::now();
    x = simpleIteration(A, b, eps, iterations);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Simple Iteration Method: " << endl;
    printVector(x);
    cout << "Iterations: " << iterations << ", Time: " << duration.count() << " microseconds" << endl;

    // Jacobi Method
    start = high_resolution_clock::now();
    x = jacobiMethod(A, b, eps, iterations);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Jacobi Method: " << endl;
    printVector(x);
    cout << "Iterations: " << iterations << ", Time: " << duration.count() << " microseconds" << endl;

    // Gauss-Seidel Method
    start = high_resolution_clock::now();
    x = gaussSeidelMethod(A, b, eps, iterations);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Gauss-Seidel Method: " << endl;
    printVector(x);
    cout << "Iterations: " << iterations << ", Time: " << duration.count() << " microseconds" << endl;

    return 0;
}