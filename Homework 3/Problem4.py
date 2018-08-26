from __future__ import division
from __future__ import print_function
from pylab import *


def trapezoidal(f, a, b, N):
    area = 0
    for i in range(N):
        x_left = a + i * (b-a)/N
        x_right = a + (i+1) * (b-a)/N
        height = (b-a)/N
        base1 = f(x_right)
        base2 = f(x_left)
        area += (base1 + base2)*height/2
    return area


def midpoint(f, a, b, N):
    area = 0
    for i in range(N):
        x_left = a + i * (b-a)/N
        x_right = a + (i+1) * (b-a)/N
        diff = x_right - x_left
        area += diff * f(x_left + diff/2)
    return area

if __name__ == "__main__":

    # Variable setup

    correct_value = 2
    left_bound = 0
    right_bound = pi
    function = sin

    N_max = 100
    N_min = 2

    h_list = zeros(N_max - N_min + 1)
    errors = array([zeros(N_max - N_min + 1), zeros(N_max - N_min + 1)])

    for N in range(N_min, N_max+1):
        print("Calculating for N =", N)
        h = (right_bound - left_bound)/N
        h_list[N-N_min] = h
        errors[0][N-N_min] = fabs(correct_value - trapezoidal(function, left_bound, right_bound, N))
        errors[1][N-N_min] = fabs(correct_value - midpoint(function, left_bound, right_bound, N))
        print("Error for trapezoidal:", errors[0][N-N_min])
        print("Error for midpoint:", errors[1][N-N_min], "\n")

    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.gca().invert_xaxis()

    plot(h_list, errors[0], 'b-',
         h_list, errors[1], 'g-',
         h_list, h_list * h_list, 'r-')
    xlabel("h value", fontsize=20)
    ylabel("error", fontsize=20)
    legend(["Trapezoidal", "Midpoint", "$h^2$"],
           loc="upper center")
    title("error vs h value", fontsize=20)

    plt.show()
