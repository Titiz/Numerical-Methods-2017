from __future__ import division
from __future__ import print_function
from pylab import *


def f(x):
    return cos(4*sin(x/2)*pi)


def f_prime(x):
    return -2 * pi * sin(4 * pi * sin(x/2)) * cos(x/2)


def derivative(f, points):
    f_values = zeros(points)
    fp_values = zeros(points)
    for k in range(points):
        current_angle = 2 * pi / points * k
        f_values[k] = f(current_angle)
        fp_values[k] = f_prime(current_angle)
    f_values = rfft(f_values)
    for k in range(points//2 + 1):
        f_values[k] *= (k*1j)
    return irfft(f_values)


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
        for j in range(point_number):
            abs_value = abs(correct_values[i][j] - derivative_values[i][j])
            if max_abs_value < abs_value:
                max_abs_value = abs_value

        abs_values[i] = max_abs_value


    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')

    for i in range(len(point_numbers)):
        plot(point_numbers, abs_values, "r")

    xlabel("number of discrete points", fontsize=20)
    ylabel("$error$", fontsize=20)

    plt.show()


