# I'm using Python 2.x.
# The following import makes division work as in Python 3. Thus 2/3=0.66666....
# rather than 2/3=0 (in Python 2 division between integers yields an integer;
# Python 3 has the additional operator // for integer division)
from __future__ import division
from __future__ import print_function
from pylab import *
from scipy.linalg import solve


def PolyLagrange(absc, ordi, x):
    """Evaluates at x the polynomial passing through the points having
abscissas 'absc' and ordinates 'ordi' using Lagrange's method.

Input:

absc, ordi: 1D arrays or iterable python objects that convert to 1D
            array. They must have the same length.

x: 1D array or a number

Output: 1D array of same length as x, or a number if x is a number.

    """
    absc = array(absc)
    ordi = array(ordi)
    if shape(absc) != shape(ordi):
        raise ValueError("Abscissas and ordinates have different length!")
    if len(shape(absc)) != 1:
        raise ValueError("Abscissas and ordinates must be 1D objects!")
    Npts = len(absc)
    product = ones((Npts,) + shape(x))
    for i in range(Npts):
        for j in range(Npts):
            if i == j: continue
            product[i] *= (x - absc[j]) / (absc[i] - absc[j])
    return dot(ordi, product)


# ------------------------------------------------------------
# Hints for exercise n.3
def SplineCoefficients(absc, ordi):
    """Computes a matric containing the coefficients of the polynomials
forming a set of natural cubic splines interpolating the points (absc, ordi).

Input:

absc, ordi: 1D arrays or iterable python objects that convert to 1D
            array. They must have the same length. absc must be ordered
            in growing order.

Output: 

A matrix with four columns and as many rows as the elements in absc
minus one. If c0, c1, c2, c3 are the elements along the i-th row, the
corresponding interpolating polynomial is 

S_i(x)=c0 + c1*(x - absc[i]) + c2*(x - absc[i])**2 + c3*(x - absc[i])**3

    """
    eq_size = (len(absc))
    equation_matrix = zeros((eq_size, eq_size))  # Matrix full of 0s is created
    equation_matrix[0][0] = 1  # We set the first term of the diagonal = 1
    equation_matrix[-1][-1] = 1 # We set the last term of the diagonal = 1
    right_matrix = [0] * eq_size  # We create the matrix on the right side of the equation
    for i in range(eq_size-2):  # We loop through the rest of the rows
        delta_x = absc[i+1] - absc[i]
        delta_y = ordi[i+1] - ordi[i]
        delta_x_next = absc[i+2] - absc[i+1]
        delta_y_next = ordi[i+2] - ordi[i+1]

        equation_matrix[i+1][i] = delta_x  # We change three entries in every row
        equation_matrix[i+1][i+1] = 2*(delta_x + delta_x_next)
        equation_matrix[i+1][i+2] = delta_x_next
        # We also update the matrix that is on the right side of the equality.
        right_matrix[i+1] = 3*(delta_y_next /delta_x_next - delta_y / delta_x)

    result_matrix = solve(matrix(equation_matrix), matrix(right_matrix).transpose())  # We solve the system

    coefficient_matrix = zeros((eq_size-1, 4))  # We then create our 4 x (n-1) matrix full of zeros

    for i in range(eq_size-1):  # We now go through every row and update the values of the coefficients
        delta_x = absc[i+1] - absc[i]
        delta_y = ordi[i+1] - ordi[i]

        coefficient_matrix[i][0] = ordi[i]
        coefficient_matrix[i][1] = delta_y/delta_x - (delta_x)/3 * (2*result_matrix[i] + result_matrix[i+1])
        coefficient_matrix[i][2] = result_matrix[i]
        coefficient_matrix[i][3] = (result_matrix[i+1] - result_matrix[i])/(3 * (absc[i+1] - absc[i]))

    return coefficient_matrix  # We finally return the matrix that will be used for evaluating splines



def __evaluate_spline(coeff, absc, x):
    """Don't use this directly. Call 'evaluate_spline', instead."""

    if x <= absc[0]:  # left interpolation
            index = 1

    elif x >= absc[-1]:  # right interpolation
        index = len(absc)-1

    else:   # the x falls in the range of our further left and right abscissas
        index = searchsorted(absc, x)

    result = 0  # where the y value of the of the polynomial at x will be stored

    for i in range(len(coeff[index-1])):
        result += coeff[index-1][i] * pow(x-absc[index-1], i)  # Calculation of the y value

    return result


def evaluate_spline(coeff, absc, x):
    """Evaluates at x the the natural spline interpolant.

Input:

coeff:  matrix containinig the coefficients of the cubic interpolants.
        This is the output of 'SplineCoefficients'.

absc:   the same abscissas passed to 'SplineCoefficients' to compute 'coeff'.

x:      1D array or a number

Output: 1D array of same length as x, or a number if x is a number.

    """
    # x is a number
    if shape(x) == ():
        p = __evaluate_spline(coeff, absc, x)
    # x is a 1D array, or list, or tuple
    else:
        x = array(x)
        p = zeros_like(x)
        for i in xrange(len(x)):
            p[i] = __evaluate_spline(coeff, absc, x[i])
    return p


def f(x):
    return 1/(1+x*x)


# ------------------------------------------------------------

if __name__ == "__main__":

    N = 9  # number of points
    interval = [-5, 5]
    spacing = 10 / N

    points = [interval[0] + spacing / 2 + spacing * m for m in range(N)]
    values = [f(x) for x in points]

    print("The interpolation points are as follows:")
    print(points)

    print("Respectively, the function value at these points are:")
    print(values)

    coefficients = SplineCoefficients(points, values)

    max_diff = 0
    max_x = 0
    poly_values = zeros(1001)
    x = linspace(interval[0], interval[1], 1001)

    for i in range(1001):
        x_coord = x[i]
        y_coord = evaluate_spline(coefficients, points, x_coord)
        fx = f(x_coord)
        if abs(fx - y_coord) > max_diff:
            max_diff = abs(fx - y_coord)
            max_x = x_coord
        poly_values[i] = y_coord

    print("The maximum error found on the test points was:", max_diff)
    print("This maximum error was found at the x value:", max_x)

    print("Plotting the graph...")

    plot(x, poly_values, 'b-',
         x, f(x), 'g-',
         points, values, 'or')
    xlabel("Abscissas", fontsize=20)
    ylabel("Ordinates", fontsize=20)
    legend(["Polynomial fit", "Original function", "Interpolation nodes"],
           loc="upper center")
    title("Maximum absolute error: {:7.4f}".format(max_diff), fontsize=20)

    plt.show()
