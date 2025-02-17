#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <windows.h>

// Установка кодировки консоли на UTF-8
void setConsoleToUTF8() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
}

// Класс плотной матрицы
class DenseMatrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows, cols;

public:
    DenseMatrix(const std::vector<std::vector<double>>& input) : data(input) {
        rows = data.size();
        cols = (rows > 0) ? data[0].size() : 0;
        for (const auto& row : data) {
            if (row.size() != cols) {
                throw std::invalid_argument("All rows must have the same length.");
            }
        }
    }

    double operator()(size_t i, size_t j) const {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("Index out of range.");
        }
        return data[i][j];
    }

    std::vector<double> operator*(const std::vector<double>& vec) const {
        if (vec.size() != cols) {
            throw std::invalid_argument("Vector size must match the number of columns.");
        }
        std::vector<double> result(rows, 0.0);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result[i] += data[i][j] * vec[j];
            }
        }
        return result;
    }

    void print() const {
        for (const auto& row : data) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << "\n";
        }
    }
};

// Класс CSR-матрицы
class CSRMatrix {
private:
    std::vector<double> values;
    std::vector<size_t> col_indices;
    std::vector<size_t> row_ptr;
    size_t rows, cols;

public:
    CSRMatrix(const std::map<std::pair<size_t, size_t>, double>& dok, size_t num_rows, size_t num_cols)
        : rows(num_rows), cols(num_cols) {
        row_ptr.push_back(0);
        size_t current_row = 0;

        for (const auto& entry : dok) {
            size_t i = entry.first.first;
            size_t j = entry.first.second;
            double val = entry.second;

            while (current_row < i) {
                row_ptr.push_back(values.size());
                current_row++;
            }

            values.push_back(val);
            col_indices.push_back(j);
        }

        while (row_ptr.size() <= rows) {
            row_ptr.push_back(values.size());
        }
    }

    std::vector<double> operator*(const std::vector<double>& vec) const {
        if (vec.size() != cols) {
            throw std::invalid_argument("Vector size must match the number of columns.");
        }
        std::vector<double> result(rows, 0.0);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
                result[i] += values[j] * vec[col_indices[j]];
            }
        }
        return result;
    }

    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
                std::cout << "(" << i << ", " << col_indices[j] << ") = " << values[j] << "\n";
            }
        }
    }
};

// Генерация случайной плотной матрицы
DenseMatrix generate_random_dense_matrix(size_t rows, size_t cols, double sparsity = 0.0) {
    std::vector<std::vector<double>> data(rows, std::vector<double>(cols, 0.0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (dist(gen) > sparsity) {
                data[i][j] = dist(gen);
            }
        }
    }
    return DenseMatrix(data);
}

// Генерация случайной CSR-матрицы
CSRMatrix generate_random_csr_matrix(size_t rows, size_t cols, double sparsity = 0.0) {
    std::map<std::pair<size_t, size_t>, double> dok;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (dist(gen) > sparsity) {
                dok[{i, j}] = dist(gen);
            }
        }
    }
    return CSRMatrix(dok, rows, cols);
}

// Измерение времени выполнения
template<typename Func>
double measure_time(Func func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

// Экспорт результатов в файл
void export_results(const std::vector<size_t>& sizes, const std::vector<double>& dense_times, const std::vector<double>& csr_times) {
    std::ofstream file("results.csv");
    file << "Razmer,Plotnaya matrica,CSR matrica\n"; // Заголовки столбцов
    for (size_t i = 0; i < sizes.size(); ++i) {
        file << sizes[i] << "," << dense_times[i] << "," << csr_times[i] << "\n"; // Разделитель - запятая
    }
    file.close();
}

int main() {
    setConsoleToUTF8(); // Устанавливаем кодировку консоли

    std::vector<size_t> sizes = { 100, 500, 1000, 2000 };
    std::vector<double> dense_times, csr_times;

    for (size_t size : sizes) {
        std::vector<double> vec(size, 1.0);

        auto dense_matrix = generate_random_dense_matrix(size, size, 0.9);
        auto csr_matrix = generate_random_csr_matrix(size, size, 0.9);

        double dense_time = measure_time([&]() { auto result = dense_matrix * vec; });
        double csr_time = measure_time([&]() { auto result = csr_matrix * vec; });

        dense_times.push_back(dense_time);
        csr_times.push_back(csr_time);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Razmer matricy: " << size << "x" << size << "\n";
        std::cout << "Vremya umnozheniya plotnoy matricy: " << dense_time << " sec.\n";
        std::cout << "Vremya umnozheniya CSR-matricy: " << csr_time << " sec.\n";
        std::cout << "----------------------------------------\n";
    }

    export_results(sizes, dense_times, csr_times);

    return 0;
}