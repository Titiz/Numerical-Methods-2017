from __future__ import division
from __future__ import print_function
from math import pi
from pylab import *

import Problem1 as pb1

def chebyshev(a, b, N):
    """Caculates N Chebyshev zeros on the interval [a,b]

    Input:

    a:      left bound of the interval

    b:      right bound of the interval

    N:      number of Chebyshev zeros wanted to be found

    Output:

    zeros:  list containing N Chebyshev zeros

        """
    zeros = []  # This is the list we will return
    alpha = (b+a)/2  # Quantities for changing bounds
    beta  = (b-a)/2
    for i in range(1, N+1):  # We calculate N chebyshev zeros
        zero = alpha + beta*cos(pi*(2*i-1)/(2*N))
        zeros.append(zero)
    return zeros


def f(x):
    return 1/(1+x*x)

if __name__ == "__main__":

    N = 9  # number of points
    interval = [-5, 5]
    cheb_zeros = chebyshev(interval[0], interval[1], N)  # retrieve the chebyshev zeros
    values = [f(x) for x in cheb_zeros]  # calculate the function for the zeros

    print("The interpolation points are as follows:")
    print(cheb_zeros)

    print("Respectively, the function value at these points are:")
    print(values)

    P = pb1.newton_interpolation(values, cheb_zeros)
    X = cheb_zeros

    print("This gives the following interpolation polynomial:")
    pb1.print_polynomial(P, X)

    max_diff = 0
    max_x = 0
    poly_values = zeros(1001)
    x = linspace(interval[0], interval[1], 1001)

    for i in range(1001):
        x_coord = x[i]
        y_coord = pb1.get_value_of_polynomial(P, X, x_coord)
        fx = f(x_coord)
        if abs(fx - y_coord) > max_diff:
            max_diff = abs(fx - y_coord)
            max_x = x_coord
        poly_values[i] = y_coord

    print("The maximum error found on the test points was:", max_diff) # max difference
    print("This maximum error was found at the x value:", max_x)

    print("Plotting the graph...")

    figure(1)
    plot(x, poly_values, 'b-',
         x, f(x), 'g-',
         cheb_zeros, values, 'or')
    xlabel("Abscissas", fontsize=20)
    ylabel("Ordinates", fontsize=20)
    legend(["Polynomial fit", "Original function", "Interpolation nodes"],
           loc="upper center")
    title("Maximum absolute error: {:7.4f}".format(max_diff), fontsize=20)

    plt.show()
