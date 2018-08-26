from __future__ import print_function, division
from pylab import *
from main import *
''' This file is used for graphing the error with respect to time steps using coeffs instead of expanded 0s.
Instead of creating a polynomial out of zeros we create a polynomial out of coefficients and
then we check if our calculated zeros evaluate to 0 when plugged into the polynomial'''
if __name__ == "__main__":
    lowest_order = 2 # Lowest order of tested polynomial, I always use 2.
    highest_order = 15 # Highest order of test polynomial
    coeff_bound = 10000 # Maximum coefficients of the randomly generated polynomial
    pol_bound = 10000 # Maximum initial guess`
    tests_per_order = 10 # How many tests should be done per 1 order of a polynomial
    h = 0.01  # what is the time-step
    iteration_numbers = [10, 15, 20]  # what iteration numbers should be used.
    T = 1 # T-value described in the algorithm. Always 1 for my purposes.

    for iteration_number in iteration_numbers:
        average_errors = zeros(highest_order - lowest_order + 1)  # array storing largest errors, for plotting purposes
        degrees = array(range(lowest_order, highest_order + 1))  # array storing degree, for plotting purposes.
        for degree in range(lowest_order, highest_order+1):  # We iterate through each degree.
            largest_error = 0  # Will store the largest error of the polynomial
            total_error = 0 # Will store the cumulative errors of all trials. Later used to find mean error.
            for j in range(tests_per_order):
                print("i : j |", degree, ':', j) # Printed so user can see progress of the program

                # We generate random coefficients from which we make a polynomial
                coeffs = coeff_bound*(rand(degree) + 1j * rand(degree))
                p = Polynomial(degree, coeffs, [-pol_bound, pol_bound], False)

                # We perform the iterative numerical integration.
                for i in range(iteration_number):
                    p.evaluate_time_zeros(h, T)

                # We look at the error in  the coefficients
                approx_coeffs = generate_coeffs(degree, p.time_zeros)
                difference = coeffs - approx_coeffs
                error = norm(difference)

                # We check if the trial had an error larger than any previous trial
                if error > largest_error:
                    largest_error = error # if so, we change the largest error
                # Regardless, we add the biggest error of this trial to the total error.
                total_error += error

            average_errors[degree - lowest_order] = total_error / tests_per_order # Calculate evarge error for degree

        semilogy(degrees, average_errors)
    title('iteration comparison')
    xlabel("polynomial degree")
    ylabel("average error")
    legend(["10", "15", "20"], loc='lower left')
    show()









