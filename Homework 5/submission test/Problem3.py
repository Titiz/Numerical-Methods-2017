from __future__ import print_function
from __future__ import division
from pylab import *
from scipy.linalg import *
from Problem1 import create_matrix
from Problem4 import inv_lower_triangular # This is a function that returns the inverse of a lower triangular matrix


def rel_back_error(A, x, b):
    return norm(b - dot(A, x), ord=inf)/norm(b, ord=inf)


def gauss_schneidel(A, b, x=None, iterations=30, threshold=1e-12):
    if x is None:
        x = ones(len(A))

    D = diag(diag(A))
    U = triu(A, 1)
    L = tril(A, -1)
    DL_inv = inv_lower_triangular(D, L)

    iteration = 0
    while iteration < iterations and rel_back_error(A, x, b) > threshold:
        x = dot(DL_inv, b - dot(U, x))
        iteration += 1

    return x


if __name__ == "__main__":

    print("Calculating...")

    n = 100

    matrix = create_matrix(n)

    b = rand(n)

    x = gauss_schneidel(matrix, b)

    print("The solution vector is:\n", x)

    print("The relative backwards error is:", rel_back_error(matrix, x, b))
