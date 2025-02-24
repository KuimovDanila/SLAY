#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Function to compute the norm of a vector
double norm(const vector<double>& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x * x;
    }
    return sqrt(sum);
}

// Function to compute the sign of a number
double sign(double x) {
    return (x >= 0) ? 1.0 : -1.0;
}

// Function to compute the Householder vector
vector<double> householder(vector<double>& x) {
    int n = x.size();
    vector<double> v(n, 0.0);
    double sigma = 0.0;

    for (int i = 1; i < n; ++i) {
        sigma += x[i] * x[i];
    }

    if (sigma == 0.0) {
        v[0] = 1.0;
        return v;
    }

    double mu = sqrt(x[0] * x[0] + sigma);
    double beta = x[0] + sign(x[0]) * mu;

    for (int i = 1; i < n; ++i) {
        v[i] = x[i] / beta;
    }

    v[0] = 1.0;
    return v;
}

// Function to apply Householder transformation to a matrix
void apply_householder(vector<vector<double>>& A, const vector<double>& v, int k) {
    int n = A.size();
    int m = A[0].size();
    double beta = -2.0 / (norm(v) * norm(v));

    for (int j = k; j < m; ++j) {
        double sum = 0.0;
        for (int i = k; i < n; ++i) {
            sum += v[i - k] * A[i][j];
        }
        sum *= beta;
        for (int i = k; i < n; ++i) {
            A[i][j] += sum * v[i - k];
        }
    }
}

// Function to perform QR decomposition
void qr_decomposition(vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size();
    int m = A[0].size();
    Q = vector<vector<double>>(n, vector<double>(n, 0.0));
    R = A;

    // Initialize Q as identity matrix
    for (int i = 0; i < n; ++i) {
        Q[i][i] = 1.0;
    }

    for (int k = 0; k < m; ++k) {
        vector<double> x(n - k);
        for (int i = k; i < n; ++i) {
            x[i - k] = R[i][k];
        }

        vector<double> v = householder(x);
        apply_householder(R, v, k);

        // Apply Householder transformation to Q
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int i = k; i < n; ++i) {
                sum += v[i - k] * Q[i][j];
            }
            sum *= -2.0 / (norm(v) * norm(v));
            for (int i = k; i < n; ++i) {
                Q[i][j] += sum * v[i - k];
            }
        }
    }
}

// Function to solve a linear system using QR decomposition
vector<double> solve_qr(const vector<vector<double>>& Q, const vector<vector<double>>& R, const vector<double>& b) {
    int n = Q.size();
    vector<double> y(n, 0.0);
    vector<double> x(n, 0.0);

    // Compute y = Q^T * b
    for (int i = 0; i < n; ++i) {
        y[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            y[i] += Q[j][i] * b[j];
        }
    }

    // Back substitution to solve Rx = y
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}

int main() {
    // Example matrix A and vector b
    vector<vector<double>> A = {
        {12, -51, 4},
        {6, 167, -68},
        {-4, 24, -41}
    };

    vector<double> b = { 1, 2, 3 };

    // Perform QR decomposition
    vector<vector<double>> Q, R;
    qr_decomposition(A, Q, R);

    // Solve the system Ax = b
    vector<double> x = solve_qr(Q, R, b);

    // Output the results
    cout << "Matrix Q:" << endl;
    for (const auto& row : Q) {
        for (double x : row) {
            cout << x << " ";
        }
        cout << endl;
    }

    cout << "Matrix R:" << endl;
    for (const auto& row : R) {
        for (double x : row) {
            cout << x << " ";
        }
        cout << endl;
    }

    cout << "Solution x:" << endl;
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << endl;

    return 0;
}