from __future__ import division
from __future__ import print_function
from pylab import *


def derivative_1(f, x, h):
    return (f(x+h) - f(x))/h


def derivative_2(f, x, h):
    return (f(x+h) - f(x-h))/(2*h)


def derivative_3(f, x, h):
    return (-2*f(x-h) - 3*f(x) + 6*f(x+h) - f(x+2*h))/(6*h)


def derivative_4(f, x, h):
    return (f(x-2*h) - 8 * f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h)


def f(x):
    return exp(x)


if __name__ == "__main__":

    iteration_number = 75
    correct_value = 1  # Derivative exp at 0

    print("Calculating...")

    max_errors = array([zeros(iteration_number),
                        zeros(iteration_number),
                        zeros(iteration_number),
                        zeros(iteration_number)])

    h_list = zeros(iteration_number)

    for i in range(iteration_number):
        h = pow(10, -i/5)
        h_list[i] = h
        print("Iteration =", i, "h =", h)
        derivatives = array([derivative_1(f, 0, h), derivative_2(f, 0, h), derivative_3(f, 0, h), derivative_4(f, 0, h)])
        for j in range(len(derivatives)):
            max_errors[j][i] = fabs(correct_value - derivatives[j])
            print("Error for derivative_" + str(j)+":", fabs(derivatives[j] - correct_value))

    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.gca().invert_xaxis()

    plot(h_list, max_errors[0], 'b-',
         h_list, max_errors[1], 'g-',
         h_list, max_errors[2], 'r-',
         h_list, max_errors[3], 'c-')
    xlabel("h value", fontsize=20)
    ylabel("error", fontsize=20)
    legend(["$P_1$", "$P_2$", "$P_3$", "$P_4$"],
           loc="upper center")
    title("error vs h value", fontsize=20)

    plt.show()