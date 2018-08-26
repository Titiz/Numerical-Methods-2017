from __future__ import print_function
from __future__ import division
from pylab import *
from scipy.linalg import *
from Problem1 import create_matrix


def rel_back_error(A, x, b):
    return norm(b - dot(A, x), ord=inf)/norm(b, ord=inf)


def jacobi_iteration(A, b, x, threshold = 1e-12):
    D = diag(A)

    D_inv = array([1/D[k] for k in range(len(D))])
    D_inv = diag(D_inv)

    R = A - diag(D)

    while rel_back_error(A, x, b) > threshold:
        x = dot(D_inv, b) - dot(dot(D_inv, R), x)

    return x


if __name__ == "__main__":
    n = randint(3, 100)

    print("matrix size chosen to be n =", n)

    matrix1 = create_matrix(n)

    b_1 = rand(n)

    initial_guess = ones(n)

    x_1 = jacobi_iteration(matrix1, b_1, initial_guess)

    error_1 = norm(dot(matrix1, x_1) - b_1, ord=inf)/norm(b_1, ord=inf)

    print("The solution is given by the vector: \n", x_1)
    print("relative backward error for matrix1 is:", error_1)


