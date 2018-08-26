from __future__ import print_function
from __future__ import division
from pylab import *
from scipy.linalg import *


def create_matrix(n):
    """Returns an nxn diagonally-dominant matrix to be used for the
exercises.

Input:
    n  a positive integer

Output:
    an nxn matrix

Note: in the PLU decomposition of the output matrix, P is the identity matrix.
"""
    n = int(n)
    if n<3:
        raise(ValueError("n must be at least 3."))
    A = empty((n,n))#3.*identity(n)
    v = zeros((n,))
    v[0]  = 3.
    v[1]  = 1.
    v[-1] = 1.
    for i in xrange(n):
        A[i] = roll(v, i)

    return A


def create_matrix_2(n):
    """Returns an nxn matrix to be used for the exercises.

Input:
    n  a positive integer

Output:
    an nxn matrix

Note: in the PLU decomposition of the output matrix, P is not the identity matrix. The matrix returned by this function is not, in general, diagonally dominant.
"""
    A = create_matrix(n)
    primes = [3, 5, 7, 11]
    lp = len(primes)
    # swap the rows in a non-trivial way
    for i in xrange(n):
        ip = i % lp
        j = (i + primes[ip]) % n
        line = copy(A[i])
        A[i] = A[j]
        A[j] = line

    return A

# ------------------------------------------------------------------------- #

def back_substitute(LUP,  b):
    P = transpose(LUP[0])
    L = LUP[1]
    U = LUP[2]
    n = len(P)

    Pb = dot(P, b)
    c_array = zeros(n)

    for i in range(0, n):
        c_array[i] = Pb[i] - dot(L[i][:i], c_array[:i])

    x_array = zeros(n)
    for i in range(n-1, -1, -1):
        if i == 0:  # Takes care of a splicing issue, where specifying [n-1:-1:-1] does not include the 0th element.
            x_array[i] = (c_array[i] - dot(U[i][n-1::-1], x_array[n - 1:: -1])) / U[i][i]
        else:
            x_array[i] = (c_array[i] - dot(U[i][n-1:i-1:-1], x_array[n-1: i-1: -1])) / U[i][i]
    return x_array


if __name__ == "__main__":

    n = randint(3, 100)

    print("matrix size chosen to be n =", n)

    matrix1 = create_matrix(n)
    matrix2 = create_matrix_2(n)

    b_1 = rand(n)
    b_2 = rand(n)

    LUP1 = lu(matrix1)
    LUP2 = lu(matrix2)

    x_1 = back_substitute(LUP1, b_1)
    x_2 = back_substitute(LUP2, b_2)

    error_1 = norm(dot(matrix1, x_1) - b_1, ord=inf)
    error_2 = norm(dot(matrix2, x_2) - b_2, ord=inf)

    print("backwards error for matrix1 is:", error_1)
    print("backwards error for matrix2 is:", error_2)


