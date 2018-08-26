from __future__ import print_function, division
from main import *
from pylab import *
''' This file is used for graphing the error with respect to different timesteps. Here zeros were generated, expanded
and the expansion gave a polynomial with known zeros. The polynomial are solved numerically and the approximated
zeros are compared to the initially generated ones.
'''
if __name__ == "__main__":
    lowest_order = 2 # Lowest order of tested polynomial, I always use 2.
    highest_order = 15 # Highest order of test polynomial
    zero_bound = 2 # Maximum zero of the randomly generated polynomial
    pol_bound = 10 # Maximum initial guess
    tests_per_order = 10 # How many tests should be done per 1 order of a polynomial
    iteration_number = 10  # How many iterations the algorithm should perform
    hs = [0.01, 0.001] # what timesteps should be used
    T = 1 # T-value described in the algorithm. Always 1 for my purposes.

    for h in hs:
        average_errors = zeros(highest_order - lowest_order + 1)  # array storing largest errors, for plotting purposes
        degrees = array(range(lowest_order, highest_order + 1))  # array storing degree, for plotting purposes.
        for degree in range(lowest_order, highest_order+1):  # We iterate through each degree.
            largest_error = 0  # Will store the largest error of the polynomial
            total_error = 0 # Will store the cumulative errors of all trials. Later used to find mean error.
            for j in range(tests_per_order):
                print("i : j |", degree, ':', j) # Printed so user can see progress of the program

                # We generate random zeros from which we make a polynomial
                correct_zeros = zero_bound*(rand(degree) + 1j*rand(degree))
                p = Polynomial(degree, correct_zeros, [-pol_bound, pol_bound])

                # We perform the iterative numerical integration.
                for i in range(iteration_number):
                    p.evaluate_time_zeros(h, T)

                # We sort the approximated zeros and the original zeros by their real part.
                correct_zeros = sort(correct_zeros)  # sort the original zeros
                p.time_zeros = sort(p.time_zeros)  # sort the approximated zeros
                # estimate zero error by comparison to the originally generated zeros.
                difference = absolute(correct_zeros - p.time_zeros) # estimate zero error
                error = norm(difference)

                # We check if the trial had an error larger than any previous trial
                if error > largest_error:
                    largest_error = error # if so, we change the largest error
                # Regardless, we add the biggest error of this trial to the total error.
                total_error += error

            average_errors[degree - lowest_order] = total_error / tests_per_order # Calculate evarge error for degree

        semilogy(degrees, average_errors)
    title('time-step comparison')
    xlabel("polynomial degree")
    ylabel("average error")
    legend(["h = 0.01", "h= 0.001"], loc='lower right')
    show()









