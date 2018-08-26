from pylab import *
from scipy.linalg import qr

coeff = array([1./4, -4.8, 30., -70., 25.])

seed(1234567)
def noisy_quartic(c, m=100):
    """Returns a set of m points constructed by adding gaussian,
independent random values to the value of a quartic polynomial having
coefficients c evaluated at random positions in [0,10].

Input:
c:  the vector of the 5 coefficients of the quartic (c[0] being the
    highest order one)
m:  the number of points to be generated

Output:
x, y: the x, y coordinates of the m points
    """
    x = 10.*rand(m)
    y = (x*(x*(x*(x*c[0] + c[1]) + c[2]) + c[3]) + c[4] + 10*randn(m))
    return x, y

x, y = noisy_quartic(coeff)

# -- Exercise 1

m = 100
coeff_count = 5

if __name__ == "__main__":
    matrix_A = zeros((m, coeff_count))
    for i in range(m):
        for j in range(coeff_count):
            matrix_A[i, j] = pow(x[i], j)

    matrix_Q, matrix_R = qr(matrix_A)

    vector_b = y

    vector_d = dot(transpose(matrix_Q), vector_b)

    matrix_up_R = matrix_R[:coeff_count]
    vector_hat_d = vector_d[:coeff_count]

    vector_x = solve(matrix_up_R, vector_hat_d)

    figure()

    xlim([0, 10])

    plot(x,y, "o")

    t = np.arange(0, 10, 0.1)

    plot(t, (t * (t * (t * (t * coeff[0] + coeff[1]) + coeff[2]) + coeff[3]) + coeff[4]), "-g")
    plot(t, (t*(t*(t*(t*vector_x[-1] + vector_x[-2]) + vector_x[-3]) + vector_x[-4]) + vector_x[0]), "-r")

    legend(["Points ", "Original", "Least Square"], loc="upper center")

    plot()

    plt.show()



















