import matrix
import time
import math
import matplotlib.pyplot as plt

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


def save_matrix_to_file(A, filename):
    with open(filename, 'w') as file:
        for i, row in enumerate(A.values):
            file.write(' '.join(map(str, row)))
            if i != len(A.values) - 1:
                file.write('\n')


def save_vector_to_file(B, filename):
    with open(filename, 'w') as file:
        for i, value in enumerate(B.values):
            file.write(f"{value[0]:.4f}")
            if i != len(B.values) - 1:
                file.write('\n')


def A_to_latex(filename, latex_filename):
    with open(filename, 'r') as file:
        matrix_values = [list(map(float, line.split())) for line in file.readlines()]

    latex = "\\usepackage{amsmath}\n\\[\nA_{996\\times996}=\\begin{bmatrix}\n"
    for i in range(len(matrix_values)):
        if i == 7:
            latex += ("\\vdots & \\vdots & \\vdots & \\vdots & \\vdots & \\vdots & \\vdots & \\ddots & \\vdots & "
                      "\\vdots & \\vdots & \\vdots & \\vdots & \\vdots & \\vdots \\\\\n")
            continue
        if 7 < i < len(matrix_values) - 7:
            continue
        for j in range(len(matrix_values[i])):
            if j == 7:
                latex += "\\cdots "
                continue
            if j < 7:
                latex += str(int(matrix_values[i][j])) + "& "
            if j >= len(matrix_values[i]) - 7:
                latex += "& " + str(int(matrix_values[i][j]))
            if j == len(matrix_values[i]) - 1:
                latex += " \\\\\n"
    latex += "\\end{bmatrix}\n\\]"

    with open(latex_filename, 'w') as file:
        file.write(latex)


def B_to_latex(filename, latex_filename):
    with open(filename, 'r') as file:
        vector_values = [float(line) for line in file.readlines()]

    with open(latex_filename, 'w') as file:
        file.write("\\usepackage{amsmath}\n\\[\nB_{996}=\\begin{bmatrix}\n")
        for i in range(7):
            file.write(str(vector_values[i]) + " \\\\\n")
        file.write("\\vdots \\\\\n")
        for i in range(-7, 0):
            file.write(str(vector_values[i]) + " \\\\\n")
        file.write("\\end{bmatrix}\n\\]")


def forward_substitution(A, b):
    x = matrix.Matrix(b.size)
    for row in range(A.size[0]):
        sum = 0
        for col in range(row):
            sum += x.values[col][0] * A.values[row][col]
        x.values[row][0] = (b.values[row][0] - sum) / A.values[row][row]
    return x


def backward_substitution(A, b):
    x = matrix.Matrix(b.size)
    for row in range(A.size[0] - 1, -1, -1):
        sum = 0
        for col in range(A.size[0] - 1, row - 1, -1):
            sum += x.values[col][0] * A.values[row][col]
        x.values[row][0] = (b.values[row][0] - sum) / A.values[row][row]
    return x


def norm(vector):
    sum_of_squares = 0
    for row in range(vector.size[0]):
        sum_of_squares += vector.values[row][0] ** 2
    return math.sqrt(sum_of_squares)


def residuum(A, r, b):
    return A * r - b


def lu_factorization(A):
    U = A.duplicate()
    L = matrix.identity_matrix(A.size)
    for i in range(A.size[0]):
        for j in range(i + 1, A.size[0]):
            L.values[j][i] = U.values[i][j] / U.values[i][i]
            for k in range(i, A.size[0]):
                U.values[j][k] -= L.values[j][i] * U.values[i][k]
    return L, U


def jacobi(A, b):
    D = A.get_diagonal()
    H = A.get_hollow()
    r = matrix.ones_vector(b.size)
    norm_res_j = []
    second_term = forward_substitution(D, b)
    start = time.time()
    while norm(residuum(A, r, b)) > 10 ** -10:
        norm_res_j.append(norm(residuum(A, r, b)))
        first_term = forward_substitution(D * -1, H * r)
        r = first_term + second_term
    time_j = time.time() - start
    return norm_res_j, time_j


def gauss_seidel(A, b):
    D = A.get_diagonal()
    L = A.get_lower_triangular()
    U = A.get_upper_triangular()
    r = matrix.ones_vector(b.size)
    first_term = (D + L) * -1
    second_term = forward_substitution(D + L, b)
    norm_res_gs = []
    start = time.time()
    while norm(residuum(A, r, b)) > 10 ** -10:
        norm_res_gs.append(norm(residuum(A, r, b)))
        r = forward_substitution(first_term, U * r) + second_term
    time_gs = time.time() - start
    return norm_res_gs, time_gs


def lu(A, b):
    start = time.time()
    L, U = lu_factorization(A)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    norm_res_lu = norm(residuum(A, x, b))
    time_lu = time.time() - start
    return norm_res_lu, time_lu


def solve_and_print_exercise_b(method, method_name, A, B):
    residuals, time_taken = method(A, B)
    print(f"   >> Solving using {method_name} method...")
    first_residual_norm = None
    if isinstance(residuals, list):
        filtered_residuals = []
        for res in residuals:
            if res >= 10**-9:
                filtered_residuals.append(res)
            else:
                first_residual_norm = res
                break
        print(f"   >> Residual norm:", first_residual_norm)
        print(f"   >> Iterations:", len(filtered_residuals) - 1)
        if filtered_residuals and first_residual_norm <= 10 ** -9:
            print(f"      >> Method is convergent!")
        else:
            print(f"      >> Method is not convergent!")
    else:
        print(f"   >> Residual norm:", first_residual_norm)
        print(f"   >> Method is not convergent!")
    print(f"   >> Time:", time_taken)
    return residuals


def solve_and_print_exercise_c(method, method_name, A, B):
    residuals, time_taken = method(A, B)
    print(f"   >> Solving using {method_name} method...")
    first_residual_norm = None
    if isinstance(residuals, list):
        filtered_residuals = []
        for res in residuals:
            if res <= 10**9:
                filtered_residuals.append(res)
            else:
                first_residual_norm = res
                break
        print(f"   >> Residual norm:", first_residual_norm)
        print(f"   >> Iterations:", len(filtered_residuals) - 1)
        if filtered_residuals and first_residual_norm <= 10 ** -9:
            print(f"      >> Method is convergent!")
        else:
            print(f"      >> Method is not convergent!")
    else:
        print(f"   >> Residual norm:", first_residual_norm)
        print(f"   >> Method is not convergent!")
    print(f"   >> Time:", time_taken)
    return residuals

def separator(project_name):
    alternate_character = "="
    for i in range(len(project_name) + 60):
        print(alternate_character, end="")
        if alternate_character == "=":
            alternate_character = "-"
        else:
            alternate_character = "="
    print()


def plot_exercise_b_c(jacobi_residuals, gauss_seidel_residuals, filename, threshold):
    plt.rcParams['font.family'] = 'Comic Sans MS'
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(12, 5))
    plt.plot(jacobi_residuals, 'o-', label='Jacobi', color='green')
    plt.plot(gauss_seidel_residuals, 'o-', label='Gauss-Seidel', color='orange')
    plt.axhline(y=threshold, color='lightblue', linestyle='--', linewidth=4, label='Convergence threshold')
    plt.yscale('log')
    plt.xlabel('Iteration')
    plt.ylabel('Norm of Residual (log)')
    plt.title('Change in Residual Norm with Iterations')
    plt.grid(True)
    plt.legend()
    plt.savefig(filename)


def jacobi_c(A, b):
    D = A.get_diagonal()
    H = A.get_hollow()
    r = matrix.ones_vector(b.size)
    norm_res_j = []
    second_term = forward_substitution(D, b)
    start = time.time()
    while norm(residuum(A, r, b)) < 10 ** 10:
        norm_res_j.append(norm(residuum(A, r, b)))
        first_term = forward_substitution(D * -1, H * r)
        r = first_term + second_term
    time_j = time.time() - start
    return norm_res_j, time_j


def gauss_seidel_c(A, b):
    D = A.get_diagonal()
    L = A.get_lower_triangular()
    U = A.get_upper_triangular()
    r = matrix.ones_vector(b.size)
    first_term = (D + L) * -1
    second_term = forward_substitution(D + L, b)
    norm_res_gs = []
    start = time.time()
    while norm(residuum(A, r, b)) < 10 ** 10:
        norm_res_gs.append(norm(residuum(A, r, b)))
        r = forward_substitution(first_term, U * r) + second_term
    time_gs = time.time() - start
    return norm_res_gs, time_gs


def plot_exercise_e(filename, output_filename):
    with open(filename, 'r') as file:
        times = [float(line.strip()) for line in file]

    jacobi_times = times[::3]
    gauss_seidel_times = times[1::3]
    lu_times = times[2::3]

    N_values = [100, 500, 1000, 1500, 2000, 2500, 3000]

    plt.rcParams['font.family'] = 'Comic Sans MS'
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(12, 5))
    plt.plot(N_values, jacobi_times, 'o-', label='Jacobi', color='green')
    plt.plot(N_values, gauss_seidel_times, 'o-', label='Gauss-Seidel', color='orange')
    plt.plot(N_values, lu_times, 'o-', label='LU', color='lightblue')
    plt.xlabel('Matrix Size')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time of Linear Equation Calculation Methods depending on Matrix Size')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_filename)