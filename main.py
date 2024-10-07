import os

os.system('cls')

def jacobi(A, b, tol=1e-10, max_iter=100):
    n = len(A)
    x = [0.0] * n
    iter_count = 0

    for _ in range(max_iter):
        iter_count += 1
        x_new = [0.0] * n
        for i in range(n):
            s = 0.0
            for j in range(n):
                if i != j:
                    s += A[i][j] * x[j]
            x_new[i] = (b[i] - s) / A[i][i]

        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            break

        x = x_new
        
    return f"Number of iterations:{iter_count}\nSolution: {x}"

def seidel(n, a, b, tol=1e-10, max_iter=100):
    x = [0.0] * n
    iter_count = 0

    for _ in range(max_iter):
        iter_count += 1
        x_new = [0.0] * n
        for i in range(n):
            s1 = 0.0
            s2 = 0.0
            for j in range(n):
                if j < i:
                    s1 += a[i][j] * x_new[j]
                elif j > i:
                    s2 += a[i][j] * x[j]
            x_new[i] = (b[i] - s1 - s2) / a[i][i]

        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            break

        x = x_new

    return f"Number of iterations:{iter_count}\nSolution: {x}"


matrix = [
    [30.3, 0.1278, 0.1531 , 0.1784],
    [0.0975, 29.4, 0.1481, 0.1734],
    [0.0925, 0.1178, 28.5, 0.1684],
    [0.0875, 0.1128, 0.1381, 27.6]
]
vector = [80.1684, 83.5730, 86.6095, 89.2778]

solutionjacobi = jacobi(matrix, vector)
print(solutionjacobi)

solutionseidel = seidel(4, matrix, vector)
print(solutionseidel)