import matrix
import math
import functions
import os


def exercise_a():
    print(">> [INFO] Exercise A")
    index = 193696
    e = index // 100 % 10  # 6
    f = index // 1000 % 10  # 3
    a1 = 5 + e
    a2 = a3 = -1
    N = 900 + index % 100  # 996
    print("   >> Index:", index, "| e:", e, "| f:", f, "| a1:", a1, "| a2:", a2, "| a3:", a3, "| N:", N)

    matrix_values = [a1, a2, a3]
    A = matrix.get_matrix((N, N), matrix_values)

    B = matrix.Matrix((N, 1))
    for n in range(1, N + 1):
        B.values[n - 1][0] = math.sin(n * (f + 1))

    if not os.path.exists('exercise_A'):
        os.makedirs('exercise_A')

    functions.save_matrix_to_file(A, 'exercise_A/A.txt')
    functions.save_vector_to_file(B, 'exercise_A/B.txt')
    functions.A_to_latex('exercise_A/A.txt', 'exercise_A/A_latex.tex')
    functions.B_to_latex('exercise_A/B.txt', 'exercise_A/B_latex.tex')
    print(">> [INFO] Exercise A done!")
    return A, B


def exercise_b(A, B, project_name):
    print(">> [INFO] Exercise B")
    jacobi_residuals = functions.solve_and_print_exercise_b(functions.jacobi, "Jacobi", A, B)
    functions.separator(project_name)
    gauss_seidel_residuals = functions.solve_and_print_exercise_b(functions.gauss_seidel, "Gauss-Seidel", A, B)

    functions.plot_exercise_b_c(jacobi_residuals, gauss_seidel_residuals, 'exercise_B.png', 10**-9)
    print(">> [INFO] Exercise B done!")


def exercise_c(project_name):
    print(">> [INFO] Exercise C")
    index = 193696
    a1 = 3
    a2 = a3 = -1
    N = 900 + index % 100  # 996
    print("   >> Index:", index, "| a1:", a1, "| a2:", a2, "| a3:", a3, "| N:", N)

    matrix_values = [a1, a2, a3]
    A2 = matrix.get_matrix((N, N), matrix_values)

    B2 = matrix.Matrix((N, 1))
    for n in range(1, N + 1):
        B2.values[n - 1][0] = math.sin(n * (index // 1000 % 10 + 1))

    if not os.path.exists('exercise_C'):
        os.makedirs('exercise_C')

    functions.save_matrix_to_file(A2, 'exercise_C/A.txt')
    functions.save_vector_to_file(B2, 'exercise_C/B.txt')
    functions.A_to_latex('exercise_C/A.txt', 'exercise_C/A_latex.tex')
    functions.B_to_latex('exercise_C/B.txt', 'exercise_C/B_latex.tex')

    jacobi_residuals = functions.solve_and_print_exercise_c(functions.jacobi_c, "Jacobi", A2, B2)
    functions.separator(project_name)
    gauss_seidel_residuals = functions.solve_and_print_exercise_c(functions.gauss_seidel_c, "Gauss-Seidel", A2, B2)

    functions.plot_exercise_b_c(jacobi_residuals, gauss_seidel_residuals, 'exercise_C.png', 10**9)
    print(">> [INFO] Exercise C done!")
    return A2, B2


def exercise_d(A, B, A2, B2):
    print(">> [INFO] Exercise D")
    norm_res_lu, time_lu = functions.lu(A, B)
    print("Norm of residual:", norm_res_lu)
    print("Time:", time_lu)
    norm_res_lu, time_lu = functions.lu(A2, B2)
    print("Norm of residual:", norm_res_lu)
    print("Time:", time_lu)
    print(">> [INFO] Exercise D done!")


def exercise_e():
    print(">> [INFO] Exercise E")
    index = 193696
    e = index // 100 % 10  # 6
    f = index // 1000 % 10  # 3
    a1 = 5 + e
    a2 = a3 = -1
    print("   >> Index:", index, "| e:", e, "| f:", f, "| a1:", a1, "| a2:", a2, "| a3:", a3)

    sizes = [100, 500, 1000, 1500, 2000, 2500, 3000]
    file = open("exercise_E.txt", "w")
    matrix_values = [a1, a2, a3]
    for size in sizes:
        A = matrix.get_matrix((size, size), matrix_values)
        B = matrix.Matrix((size, 1))
        for n in range(1, size + 1):
            B.values[n - 1][0] = math.sin(n * (f + 1))

        _, time_j = functions.jacobi(A, B)
        _, time_gs = functions.gauss_seidel(A, B)
        _, time_lu = functions.lu(A, B)
        file.write(str(time_j) + "\n" + str(time_gs) + "\n" + str(time_lu) + "\n")
    file.close()

    functions.plot_exercise_e('exercise_E.txt', 'exercise_E.png')
    print(">> [INFO] Exercise E done!")


def main(project_name):
    #  functions.test()
    functions.separator(project_name)
    print(project_name)
    functions.separator(project_name)
    #A, B = exercise_a()
    functions.separator(project_name)
    #exercise_b(A, B, project_name)
    functions.separator(project_name)
    #A2, B2 = exercise_c(project_name)
    functions.separator(project_name)
    #exercise_d(A, B, A2, B2)
    functions.separator(project_name)
    exercise_e()


if __name__ == '__main__':
    main('>> [TITLE] Systems of linear equations solver')