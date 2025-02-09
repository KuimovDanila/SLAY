#include "TridiagonalSolver.h"
#include <stdexcept>

TridiagonalSolver::TridiagonalSolver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
    : a(a), b(b), c(c), d(d) {
    if (a.size() != b.size() || b.size() != c.size() || c.size() != d.size()) {
        throw std::invalid_argument("All input vectors must have the same size.");
    }
}

std::vector<double> TridiagonalSolver::solve() {
    size_t n = b.size();
    std::vector<double> c_prime(n, 0.0);
    std::vector<double> d_prime(n, 0.0);
    std::vector<double> x(n, 0.0);

    // Прямой ход
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (size_t i = 1; i < n; ++i) {
        double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = c[i] * m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
    }

    // Обратный ход
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}





