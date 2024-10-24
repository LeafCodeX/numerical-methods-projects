import matrix
import math

def test():
    matrix1 = matrix.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    matrix2 = matrix.Matrix([[9, 8, 7], [6, 5, 4], [3, 2, 1]])

    print("Matrix 1:")
    for row in matrix1.values:
        print(row)

    print("\nMatrix 2:")
    for row in matrix2.values:
        print(row)

    print("\nAddition:")
    result_add = matrix1 + matrix2
    for row in result_add.values:
        print(row)

    print("\nSubtraction:")
    result_sub = matrix1 - matrix2
    for row in result_sub.values:
        print(row)

    print("\nMatrix Multiplication:")
    result_mul = matrix1 * matrix2
    for row in result_mul.values:
        print(row)

    print("\nScalar Multiplication:")
    result_scalar_mul = matrix1 * 2
    for row in result_scalar_mul.values:
        print(row)

    print("\nDiagonal Matrix:")
    diagonal_matrix = matrix1.get_diagonal()
    for row in diagonal_matrix.values:
        print(row)

    print("\nHollow Matrix:")
    hollow_matrix = matrix1.get_hollow()
    for row in hollow_matrix.values:
        print(row)

    print("\nLower Triangular Matrix:")
    lower_triangular_matrix = matrix1.get_lower_triangular()
    for row in lower_triangular_matrix.values:
        print(row)

    print("\nUpper Triangular Matrix:")
    upper_triangular_matrix = matrix1.get_upper_triangular()
    for row in upper_triangular_matrix.values:
        print(row)

    print("\nBand Matrix:")
    matrix_values = [1, 2, 3]
    A = matrix.get_matrix((10, 10), matrix_values)
    for row in A.values:
        print(row)

    print("\nIdentity Matrix:")
    identity_matrix_result = matrix.identity_matrix((4, 4))
    for row in identity_matrix_result.values:
        print(row)

    print("\nOnes Vector:")
    ones_vector_result = matrix.ones_vector((3, 1))
    for row in ones_vector_result.values:
        print(row)

def forward_substitution(L, b):
    n = L.size[0]
    y = matrix.Matrix((n, 1))
    for i in range(n):
        y.values[i][0] = b.values[i][0] - sum(L.values[i][j] * y.values[j][0] for j in range(i))
    return y

def backward_substitution(U, b):
    n = U.size[0]
    x = matrix.Matrix((n, 1))
    for i in range(n - 1, -1, -1):
        x.values[i][0] = (b.values[i][0] - sum(U.values[i][j] * x.values[j][0] for j in range(i + 1, n))) / U.values[i][i]
    return x


def lu_factorization(A):
    U = A.duplicate()
    L = matrix.identity_matrix(A.size)
    for i in range(A.size[0]):
        for j in range(i + 1, A.size[0]):
            L.values[j][i] = U.values[i][j] / U.values[i][i]
            for k in range(i, A.size[0]):
                U.values[j][k] -= L.values[j][i] * U.values[i][k]
    return L, U