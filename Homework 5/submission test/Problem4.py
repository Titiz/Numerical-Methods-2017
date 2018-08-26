from __future__ import print_function
from __future__ import division
from pylab import *
from scipy.linalg import *
from Problem1 import create_matrix


def rel_back_error(A, x, b):
    return norm(b - dot(A, x), ord=inf)/norm(b, ord=inf)


def inv_lower_triangular(D, L):  # returns the inverse of a lower triangular matrix
    size = len(D)
    I = identity(size)
    answer = I

    D_inv = array([1/diag(D)[k] for k in range(size)])
    D_inv = diag(D_inv)

    S = I
    for i in range(1, size+1):
        S = (-1) * dot(S, dot(D_inv, L))
        answer += S
    return dot(answer, D_inv)


def over_relaxation(A, b, w, x=None, iterations=30, threshold=1e-12):
    if x is None:
        x = ones(len(A))

    D = diag(diag(A))
    U = triu(A, 1)
    L = tril(A, -1)
    DwL_inv = inv_lower_triangular(D, w*L)

    iteration = 0
    while iteration < iterations and rel_back_error(A, x, b) > threshold:
        x = dot(DwL_inv, w*b - dot((w * U + (w-1) * D), x))
        iteration += 1

    return x


if __name__ == "__main__":

    print("Calculating...")

    n = 100

    matrix = create_matrix(n)

    b = rand(n)

    points_to_use = 20

    parameter_values = zeros(points_to_use)
    values = zeros(points_to_use)

    for i in range(points_to_use):
        parameter_values[i] = 1 + 0.20/points_to_use * i
        x = over_relaxation(matrix, b, parameter_values[i])
        values[i] = rel_back_error(matrix, x, b)

    print("Plotting  the graph...")

    figure(1)

    plt.plot(parameter_values, values, "r-")

    xlabel("$\omega$", fontsize=20)
    ylabel("rel. backwards error", fontsize=20)

    plt.show()