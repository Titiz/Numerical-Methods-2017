from __future__ import division
from __future__ import print_function
from pylab import *


def f(n, x):
    return cos(n * sin(x/2) * pi)

if __name__ == "__main__":

    point_amount = 128
    real_valued_point_amount = point_amount//2+1

    n_values = array([2, 4, 8, 16])

    abs_values = array([zeros(real_valued_point_amount) for i in range(len(n_values))])

    wave_numbers = arange(real_valued_point_amount)

    f_values = array([zeros(point_amount) for i in range(len(n_values))])

    for i in range(point_amount):
        for j in range(len(n_values)):
            current_angle = 2 * pi / point_amount * i
            f_value = f(n_values[j], current_angle)
            f_values[j][i] = f_value

    for i in range(len(f_values)):
        abs_values[i] = absolute(rfft(f_values[i]))

    print("Plotting  the graph...")

    figure(1)

    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')

    for i in range(len(n_values)):
        plot(wave_numbers, abs_values[i])

    xlabel("Wavenumber", fontsize=20)
    ylabel("$|C_k|$", fontsize=20)
    legend(["$n=2$", "$n=4$", "$n=8$", "$n=16$"],
           loc="lower left")

    plt.show()
