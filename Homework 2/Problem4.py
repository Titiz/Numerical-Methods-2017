from __future__ import division
from __future__ import print_function
from pylab import *
import Problem1 as pb1
import Problem2 as pb2
import Problem3 as pb3


def f(x):               # Redefine function for clarity
    return 1/(1+x*x)


if __name__ == "__main__":

    # We clear three files in which we will store values errors of our approximations.
    open('Newton.txt', 'w').close()
    open('Newton_cheb.txt', 'w').close()
    open('Spline.txt', 'w').close()

    upper_bound = 40
    lower_bound = 3

    N_values = range(lower_bound, upper_bound+1)

    max_errors = [zeros(upper_bound+1 - lower_bound),
                  zeros(upper_bound+1 - lower_bound),
                  zeros(upper_bound+1 - lower_bound)]

    for N in N_values:  # Loop as required by problem
        max_error = [0, 0, 0]  # Used store the maximum error of each interpolation
        max_x = [0, 0, 0]  # Used to store the abscicas at which the maximum error occurs

        # Generating the ordinates and absciccas

        interval = [-5, 5]
        spacing = 10 / N
        points = [interval[0] + spacing / 2 + spacing * m for m in range(N)]  # creates the absiccas
        values = [f(x) for x in points]  # creates the ordinates

        # Newton approximation:

        newton_poly = pb1.newton_interpolation(values, points)  # This gives us our newton_polynomial

        # Newton with Chebyshev zeros:

        cheb_zeros = pb2.chebyshev(-5, 5, N)  # Retrieves the Chebyshev zeros
        cheb_values = [f(x) for x in cheb_zeros]  # Retrieves the ordinates for the Chebyshev zeros
        cheb_newton_poly = pb1.newton_interpolation(cheb_values, cheb_zeros)  # Generates polynomial

        spline_coeff = pb3.SplineCoefficients(points, values)  # Returns our spline interpolation coefficients

        x = linspace(interval[0], interval[1], 1001)

        for i in range(1001):  # We loop through the 1001 points
            x_coord = x[i]
            # We evaluate all of the approximations at the iteration absica and store it in an array

            y = [pb1.get_value_of_polynomial(newton_poly, points, x_coord),
                 pb1.get_value_of_polynomial(cheb_newton_poly, cheb_zeros, x_coord),
                 pb3.evaluate_spline(spline_coeff, points, x_coord)]

            for j in range(len(y)): # We check if the current error of each approximation is bigger than the last error
                y_coord = y[j]
                fx = f(x_coord)
                if abs(fx - y_coord) > max_error[j]:
                    max_error[j] = abs(fx - y_coord)  # If yes, we record the new biggest error
                    max_x[j] = x_coord  # We also record the x value at which the biggest error occurred.

        for k in range(len(max_errors)):
            max_errors[k][N - 3] = max_error[k]

        # For personal observation, I stored each of the errors in a separate file for easier inspection.

        o = open("Newton.txt", "a")
        o.write(str(N) + ": " + str(max_error[0]) + "\n")
        o.close()

        o = open("Newton_cheb.txt", "a")
        o.write(str(N) + ": " + str(max_error[1]) + "\n")
        o.close()

        o = open("Spline.txt", "a")
        o.write(str(N) + ": " + str(max_error[2]) + "\n")
        o.close()

        # I also printed out a table of values for quick comparison.

        print("Data with value N =", N, ":")
        print("Method                | max_error        | x at which max_error was found")
        print("Newton                |", max_error[0], "|", max_x[0])
        print("Newton with Chebyshev |", max_error[1], "|", max_x[1])
        print("Cubic Spline          |", max_error[2], "|", max_x[2])
        print()

    print("outputting the graph...")

    plot(N_values, max_errors[0], 'b-',
         N_values, max_errors[1], 'g-',
         N_values, max_errors[2], 'r-')

    plt.gca().set_yscale('log')
    xlabel("N", fontsize=20)
    ylabel("Error", fontsize=20)
    legend(["Newton fit", "Newton Chebyshev fit", "Cubic Spline fit"],
           loc="upper center")
    title("Errors vs Interpolation points", fontsize=20)



    plt.show()
