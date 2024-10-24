import copy


class Matrix:
    def __init__(self, argument):
        if isinstance(argument, list):
            self.size = (len(argument), len(argument[0]))
            self.values = argument
        elif isinstance(argument, tuple):
            self.size = argument
            self.values = [[0] * argument[1] for _ in range(argument[0])]
        else:
            raise ValueError(">> [ERROR] Argument must be a list or a tuple")

    def __add__(self, other):
        if self.size != other.size:
            raise ValueError(">> [ERROR] Matrix sizes must be the same for addition")
        result = [[self.values[i][j] + other.values[i][j] for j in range(self.size[1])] for i in range(self.size[0])]
        return Matrix(result)

    def __sub__(self, other):
        if self.size != other.size:
            raise ValueError(">> [ERROR] Matrix sizes must be the same for subtraction")
        result = [[self.values[i][j] - other.values[i][j] for j in range(self.size[1])] for i in range(self.size[0])]
        return Matrix(result)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            if self.size[1] != other.size[0]:
                raise ValueError(">> [ERROR] Number of columns in the first matrix must be equal to the number of rows in the second matrix for matrix multiplication!")
            result = [[sum(self.values[i][k] * other.values[k][j] for k in range(self.size[1])) for j in range(other.size[1])] for i in range(self.size[0])]
            return Matrix(result)
        elif isinstance(other, (int, float)):
            result = [[self.values[i][j] * other for j in range(self.size[1])] for i in range(self.size[0])]
            return Matrix(result)
        else:
            raise ValueError(">> [ERROR] Multiplication is defined only for another Matrix object or a scalar value!")
    def copy(self):
        return Matrix(self.values)

    def duplicate(self):
        duplicate = Matrix(self.size)
        duplicate.values = copy.deepcopy(self.values)
        return duplicate

    def get_diagonal(self):
        diagonal = Matrix(self.size)
        diagonal_values = [self.values[i][i] for i in range(self.size[0])]
        for index, value in enumerate(diagonal_values):
            diagonal.values[index][index] = value
        return diagonal

    def get_hollow(self):
        hollow = self.duplicate()
        for i in range(hollow.size[0]):
            hollow.values[i][i] = 0
        return hollow

    def get_lower_triangular(self):
        lower_triangular = Matrix(self.size)
        for row in range(self.size[0]):
            for col in range(row):
                lower_triangular.values[row][col] = self.values[row][col]
        return lower_triangular

    def get_upper_triangular(self):
        upper_triangular = Matrix(self.size)
        for row in range(self.size[0]):
            for col in range(row + 1, self.size[1]):
                upper_triangular.values[row][col] = self.values[row][col]
        return upper_triangular


def get_matrix(size, values):
    if len(values) != 3:
        raise ValueError("Incorrect number of values for band matrix")
    a1, a2, a3 = values
    result = [[0] * size[1] for _ in range(size[0])]
    for i in range(size[0]):
        for j in range(size[1]):
            if i == j:
                result[i][j] = a1
            elif abs(i - j) == 1:
                result[i][j] = a2
            elif abs(i - j) == 2:
                result[i][j] = a3
    return Matrix(result)


def identity_matrix(size):
    result = [[0] * size[1] for _ in range(size[0])]
    for i in range(min(size)):
        result[i][i] = 1
    return Matrix(result)


def ones_vector(size):
    return Matrix([[1] for _ in range(size[0])])
