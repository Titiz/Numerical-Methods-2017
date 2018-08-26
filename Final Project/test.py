from __future__ import print_function, division
from main import *
from pylab import *
''' This file is used to test different starting conditions on the iterative algorithm. Here I was using
 the approach of first creating zeros, expanding them and then checking the solved zeros against the original zeros. '''
if __name__ == "__main__":
    lowest_order = 2 # Lowest order of tested polynomial, I always use 2.
    highest_order = 15 # Highest order of test polynomial
    zero_bound = 1000 # Maximum coefficients of the zero
    pol_bound = 1000 # Maximum initial guess coefficient
    tests_per_order = 10 # How many tests should be done per 1 order of a polynomial
    iteration_number = 1 # How many iterations the algorithm should perform
    h = 0.01 # timestep
    T = 1 # T-value described in the algorithm. Always 1 for my purposes.

    # The string below is used to name the file, to track what it actually is.
    string = str(zero_bound) + ' - ' + str(pol_bound) + '_' + str(iteration_number) + str(h)[1:] + '-O-' + str(highest_order) + ".txt"

    # Clear the files or create them
    f = open("trial_data." + string, "w")  # This stores summary of error of each degree
    g = open("trial_dataQ." + string, "w") # This stores detailed data of errors each trial
    f.close()
    g.close()

    average_errors = zeros(highest_order - lowest_order + 1) # array storing largest errors, for plotting purposes
    degrees = array(range(lowest_order, highest_order+1)) # array storing degree, for plotting purposes.

    for degree in range(lowest_order, highest_order+1):  # We iterate through each degree.
        g = open("trial_dataQ." + string, "a")  # We open in append mode to be able to edit the file
        f = open("trial_data." + string, "a")

        # We create a separator for each degree so it is easier to read.
        f.write("----------------------------------------" '\n')
        f.write("Degree = " + str(degree) + '\n')
        f.write("----------------------------------------" '\n')

        g.write("----------------------------------------" '\n')
        g.write("Degree = " + str(degree) + '\n')
        g.write("----------------------------------------" '\n')

        largest_error = 0  # Will store the largest error of the polynomial
        total_error = 0 # Will store the cummulative errors of all trials. Later used to find mean error.
        for j in range(tests_per_order):
            print("i : j |", degree, ':', j) # Printed so user can see progress of the program
            f.write("\n***Trial " + str(j+1) + ":\n")

            # We generate random zeros from which we make a polynomial
            correct_zeros = zero_bound * (rand(degree) + 1j * rand(degree))
            p = Polynomial(degree, correct_zeros, [-zero_bound, zero_bound])

            # We perform the iterative numerical integration.
            for i in range(iteration_number):
                p.evaluate_time_zeros(h, T)

            # We sort the approximated zeros and the original zeros by their real part.
            correct_zeros = sort(correct_zeros)  # sort the original zeros
            p.time_zeros = sort(p.time_zeros)  # sort the approximated zeros
            difference = absolute(correct_zeros - p.time_zeros) # estimate their error
            error = norm(difference)

            # We check if the trial had an error larger than any previous trial
            if norm(difference) > largest_error:
                largest_error = error # if so, we change the largest error
            # Regardless, we add the biggest error of this trial to the total error.
            total_error += error

            # Whe then output that data into a file
            f.write("Correct zeros: \n")
            f.write(str(correct_zeros) + '\n')
            f.write("Evaluated zeros: \n")
            f.write(str(p.time_zeros) + '\n')
            f.write("Error array: \n")
            f.write(str(difference) + '\n')
            f.write("Error = " + str(error))
            f.write('\n')
        f.write('\n')

        # Here we write the summary of the errors into the g file.
        average_errors[degree - lowest_order] = total_error / tests_per_order
        g.write("largest error of all trials = " + str(largest_error) + "\n")
        g.write("average error of all trials = " + str(total_error/tests_per_order) + '\n')
        g.write('\n')

        f.close()
        g.close()

    # We plot the results
    figure()
    title('number of iterations: ' + str(iteration_number))
    xlabel("polynomial degree")
    ylabel("average error")
    semilogy(degrees, average_errors)
    show()






