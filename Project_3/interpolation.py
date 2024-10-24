from matrix import Matrix
from math import cos, pi
from functions import lu_factorization, forward_substitution, backward_substitution

def linspace(start, end, num_points):
    return [start + i*(end-start)/(num_points-1) for i in range(num_points)]

def generate_chebyshev_nodes(start, end, num_points):
    return [(start + end)/2 + (end - start)/2 * cos((2*i+1)/(2*num_points)*pi) for i in range(num_points)]

def lagrange_interpolation_for_point(x, y, point):
    n = len(x)
    total_sum = 0
    for i in range(n):
        xi, yi = x[i], y[i]
        def g(i, n):
            prod = 1
            for j in range(n):
                if i != j:
                    xj = x[j]
                    prod *= (point - xj) / (xi - xj) if xi != xj else 1
            return prod
        total_sum += yi * g(i, n)
    return total_sum

def lagrange_interpolation(x, y, points):
    return [lagrange_interpolation_for_point(x, y, point) for point in points]


def spline_interpolation(x_points, y_points, points_to_interpolate):
    n = len(x_points)
    a = y_points
    h = [x_points[i + 1] - x_points[i] for i in range(n - 1)]

    A = Matrix((n, n))
    b = Matrix((n, 1))

    for i in range(1, n - 1):
        A.values[i][i - 1] = h[i - 1]
        A.values[i][i] = 2 * (h[i - 1] + h[i])
        A.values[i][i + 1] = h[i]
        b.values[i][0] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1])

    A.values[0][0] = 1
    A.values[n - 1][n - 1] = 1
    b.values[0][0] = 0
    b.values[n - 1][0] = 0

    L, U = lu_factorization(A)
    c = backward_substitution(U, forward_substitution(L, b))

    b_coeffs = [0] * (n - 1)
    d_coeffs = [0] * (n - 1)

    for i in range(n - 1):
        b_coeffs[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c.values[i + 1][0] + 2 * c.values[i][0]) / 3
        d_coeffs[i] = (c.values[i + 1][0] - c.values[i][0]) / (3 * h[i])

    result = []
    for point in points_to_interpolate:
        for i in range(n - 1):
            if x_points[i] <= point <= x_points[i + 1]:
                dx = point - x_points[i]
                result.append(
                    a[i] + b_coeffs[i] * dx + c.values[i][0] * dx ** 2 + d_coeffs[i] * dx ** 3)
                break
        else:
            result.append(0)

    return result
