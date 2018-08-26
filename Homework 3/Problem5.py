from __future__ import division
from __future__ import print_function
from pylab import *
import sys

def get_integrand(f, a, b):

    dx_over_dt = ((b-a)/2)

    def x(t):
        return 1/2 * ((b-a)*t + b + a)

    def integrand(t):
        return f(x(t)) * dx_over_dt

    return integrand


def __Gauss2(f):  # n = 2
    coeff = array([1, 1])
    root_list = array([-sqrt(1/3), sqrt(1/3)])
    return coeff[0] * f(root_list[0]) + coeff[1] * f(root_list[1])


def __Gauss3(f):  # n = 3
    coeff = array([5 / 9, 8 / 9, 5 / 9])
    root_list = array([-sqrt(3/5), 0, sqrt(3/5)])
    return coeff[0] * f(root_list[0]) + coeff[1] * f(root_list[1]) + coeff[2] * f(root_list[2])


def __Gauss4(f):  # n = 4
    coeff = array([(90 - 5 * sqrt(30)) / 180, (90 + 5 * sqrt(30)) / 180,
             (90 + 5 * sqrt(30)) / 180, (90 - 5 * sqrt(30)) / 180])
    root_list = array([-sqrt((15+2*sqrt(30))/35), -sqrt((15-2*sqrt(30))/35),
                 sqrt((15-2*sqrt(30))/35), sqrt((15+2*sqrt(30))/35)])

    return coeff[0] * f(root_list[0]) + coeff[1] * f(root_list[1]) + \
        coeff[2] * f(root_list[2]) + coeff[3] * f(root_list[3])


def Gauss(f, a, b, N):
    areas = zeros(3)
    for i in range(N):
        x_left = a + i*(b-a)/N
        x_right = a + (i+1)*(b-a)/N
        integrand = get_integrand(f, x_left, x_right)
        areas[0] += __Gauss2(integrand)
        areas[1] += __Gauss3(integrand)
        areas[2] += __Gauss4(integrand)
    return areas


if __name__ == "__main__":

    # First we will test whether our Gaussian quadrature functions
    # are able to integrate polynomials of degree n-1 exactly
    # We compare the float values by a margin of 10e-15.
    # This is close to the maximum precision that we can get.

    polynomials = array([
                    lambda x: 1,
                    lambda x: 1 + 2*x,
                    lambda x: 1 + 2*x + 3*x*x,
                    lambda x: 1 + 2*x + 3*x*x + 4*x*x*x,
                    lambda x: 1 + 2*x + 3*x*x + 4*x*x*x + 4*x*x*x*x,
                    lambda x: 1 + 2*x + 3*x*x + 4*x*x*x + 4*x*x*x*x + 5*x*x*x*x*x,
                    lambda x: 1 + 2*x + 3*x*x + 4*x*x*x + 4*x*x*x*x + 5*x*x*x*x*x + 6*x*x*x*x*x*x,
                    lambda x: 1 + 2*x + 3*x*x + 4*x*x*x + 4*x*x*x*x + 5*x*x*x*x*x + 6*x*x*x*x*x*x + 7*x*x*x*x*x*x*x])

    correct_results = array([2.0, 2.0, 4.0, 4.0, 28/5, 28/5, 256/35, 256/35, 2864/3])

    print("Testing Gaussian quadrature accuracy for polynomials of degree 2n-1")
    print("-" * 70)
    print("Testing Gaussian quadrature with n = 2:")
    for i in range(3):  # 2n-1 = 4-1 = 3
        if fabs(__Gauss2(polynomials[i]) - correct_results[i]) > pow(10, -15):  # We cannot compare floats by equality
            print("The error between the Gaussian quadrature and the correct integral value is too big.")
            print("Exiting the program...")
            sys.exit()
    print("Test passed")

    print("Testing Gaussian quadrature with n = 3:")
    for i in range(5):  # 2n-1 = 6-1 = 5
        if fabs(__Gauss3(polynomials[i]) - correct_results[i]) > pow(10, -15):  # We cannot compare floats by equality
            print("The error between the Gaussian quadrature and the correct integral value is too big.")
            print("Exiting the program...")
            sys.exit()
    print("Test passed")

    print("Testing Gaussian quadrature with n = 4:")
    for i in range(7):  # 2n-1 = 8-1 = 7
        if fabs(__Gauss4(polynomials[i]) - correct_results[i]) > pow(10, -15):  # We cannot compare floats by equality
            print("The error between the Gaussian quadrature and the correct integral value is too big.")
            print("Exiting the program...")
            sys.exit()
    print("Test passed")

    print("-" * 70)

    print("All tests passed. Continuing to the graphing stage. \n")

    # If we pass the test, we then graph what is required by the assignment.

    N_min = 2
    N_max = 100
    h_list = zeros(N_max - N_min + 1)
    errors = array([zeros(N_max - N_min + 1), zeros(N_max - N_min + 1), zeros(N_max - N_min + 1)])
    left_bound = 0
    right_bound = pi
    function = sin
    correct_value = 2

    print("Calculating...")

    for N in range(N_min, N_max+1):
        index = N-N_min
        h_list[index] = (right_bound - left_bound)/N
        areas = Gauss(function, left_bound, right_bound, N)
        errors[0][index] = fabs(areas[0] - correct_value)
        errors[1][index] = fabs(areas[1] - correct_value)
        errors[2][index] = fabs(areas[2] - correct_value)

    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.gca().invert_xaxis()

    plot(h_list, errors[0], 'b-',
         h_list, errors[1], 'g-',
         h_list, errors[2], 'c-')
    xlabel("h value", fontsize=20)
    ylabel("error", fontsize=20)
    legend(["Gauss2", "Gauss3", "Gauss4"],
           loc="upper center")
    title("error vs h value", fontsize=20)

    plt.show()
