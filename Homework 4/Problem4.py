from __future__ import division
from __future__ import print_function
from pylab import *

from Problem3 import derivative


def f(x):
    return 1/(1+(x-pi)*(x-pi))


def f_prime(x):
    return -(2*(x-pi))/(((x-pi)*(x-pi)+1)*((x-pi)*(x-pi)+1))


if __name__ == "__main__":
    point_numbers = 2 ** arange(3, 10)

    derivative_values = array([zeros(n) for n in point_numbers])
    correct_values = array([zeros(n) for n in point_numbers])
    abs_values = zeros(len(point_numbers))

    for i in range(len(point_numbers)):
        point_number = point_numbers[i]
        derivative_values[i] = derivative(f, point_number)
        correct_values[i] = f_prime(array([2*pi*s/point_number for s in range(point_number)]))

        max_abs_value = 0
        x_of_max_error = 0
        for j in range(point_number):
            abs_value = abs(correct_values[i][j] - derivative_values[i][j])
            if max_abs_value < abs_value:
                max_abs_value = abs_value
                x_of_max_error = 2 * pi * j / point_number
        abs_values[i] = max_abs_value
        print("n =", i+3, "has error of", max_abs_value, "at x =", x_of_max_error)

    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')

    for i in range(len(point_numbers)):
        plot(point_numbers, abs_values, "r")

    xlabel("number of discrete points", fontsize=20)
    ylabel("$error$", fontsize=20)

    plt.show()