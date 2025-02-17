import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных
data = pd.read_csv("C:/Users/LEGION/Documents/Visual Studio 2022/projects/Dense_CSRmatrix/results.csv", sep=",")

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(data["Razmer"], data["Plotnaya matrica"], label="Плотная матрица", marker='o')
plt.plot(data["Razmer"], data["CSR matrica"], label="CSR матрица", marker='s')
plt.xlabel("Размер матрицы")
plt.ylabel("Время (сек)")
plt.title("Сравнение скорости умножения матриц на вектор")
plt.legend()
plt.grid(True)
plt.show()
