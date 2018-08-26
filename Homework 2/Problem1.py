from __future__ import division
from __future__ import print_function
from pylab import *


def newton_interpolation(ordi, absc):  # The required function
    coefficients = ordi  # The starting matrix of the coefficients that we will alter each iteration
    for iteration in range(1, len(absc)):  # We iterate n-1 times
        new_coefficients = list(coefficients)  # We make a copy of the coefficients
        indices = range(iteration, len(absc))  # Each iteration we affect a smaller number of coefficients
        for index in indices:  # We update the necessary coefficients
            result = (coefficients[index] - coefficients[index-1]) / (absc[index] - absc[index-iteration])
            new_coefficients[index] = result  # We might be using old values of the coefficients, so we update the copy.
        coefficients = new_coefficients  # We now set the old coefficients to the new ones
    return coefficients  # We finally return the list of the coefficients of the polynomial.


def print_polynomial(poly, absc):  # prints out the given polynomial
    string = ""
    for n in range(len(poly)):
        string += str(float(poly[n]))
        for j in range(n):
            if (absc[j]) >= 0:
                string += "(x - " + str(absc[j]) + ")"
            else:
                string += "(x + " + str(absc[j])[1:] + ")"
        if n != len(poly)-1:
            string += " + "
    print(string)


def get_value_of_polynomial(poly, absc, x):  # gets the value of the polynomial at a given x value
    value = 0
    for n in range(len(poly)):
        product = poly[n]
        for j in range(n):
            product *= (x - absc[j])
        value += product
    return value


def f(x):
    return 1/(1+x*x)


if __name__ == "__main__":

    N = 9  # number of points
    interval = [-5, 5]
    spacing = (interval[1]-interval[0])/N

    points = [interval[0] + spacing/2 + spacing*m for m in range(N)]  # equivalent to linspace
    values = [f(x) for x in points]

    print("The interpolation points are as follows:")
    print(points)

    print("Respectively, the function value at these points are:")
    print(values)

    P = newton_interpolation(values, points)
    X = points

    print("This gives the following interpolation polynomial:")
    print_polynomial(P, X)

    max_diff = 0
    max_x = 0
    poly_values = zeros(1001)
    x = linspace(interval[0], interval[1], 1001)

    for i in range(len(x)):  # Calculates max difference and gets values of polynomial.
        x_coord = x[i]
        y_coord = get_value_of_polynomial(P, X, x_coord)
        fx = f(x_coord)
        if abs(fx - y_coord) > max_diff:  # Check for new max_difference
            max_diff = abs(fx - y_coord)
            max_x = x_coord
        poly_values[i] = y_coord  # Store the calculated value of polynomial at the given point

    print("The maximum error found on the test points was:", max_diff)  # max difference
    print("This maximum error was found at the x value:", max_x)

    print("Plotting  the graph...")

    figure(1)
    plot(x, poly_values, 'b-',
         x, f(x), 'g-',
         points, values, 'or')
    xlabel("Abscissas", fontsize=20)
    ylabel("Ordinates", fontsize=20)
    legend(["Polynomial fit", "Original function", "Interpolation nodes"],
           loc="upper center")
    title("Maximum absolute error: {:7.4f}".format(max_diff), fontsize=20)

    plt.show()

