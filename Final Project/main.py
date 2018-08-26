from __future__ import print_function, division
from pylab import *


def generate_initial_guesses(N, a, b, is_complex = True):
    """Generates an array of random numbers, bounded by a and b, with N entries.

    Input:

    N: integer of number of entries. For our purposes corresponds to degree of polynomial.

    a, b: floats denoting the bounds of the guess. a is the left bound and b is the right bound.

    is_complex: Boolean which tells the program whether real or complex guesses should be generated.

    Output: 1D array with N random entries bounded by a and b.

    """
    if is_complex:
        return array([(a + (b-a) * rand()) + (a + (b-a) * rand())*1j for _ in range(N)])
    else:
        return array([(a + (b-a) * rand()) for _ in range(N)])


def get_next_index_tracker(N, index_tracker):
    """ Helper function for genereate_coeffs.

     Input:

     N: size of the array

     index_tracker: an array of indices of the coefficient vector c described in the paper. It tracks all possible
     combinations of n_1, n_2, ... n_m values seen in under the sum term in (4c) of the original paper.

     Output: an updated version of the index tracker.
     """

    pointer = 0
    while True:
        if index_tracker[pointer] < N - pointer:
            index_tracker[pointer] += 1
            break
        pointer += 1
    for i in reversed(range(0, pointer)):
        index_tracker[i] = index_tracker[i+1] + 1
    return index_tracker


def generate_coeffs(N, zero_array):
    """ Take an array of zeros a polynomial and returns an array of coefficients of that polynomial.

    Input:

    N: degree of polynomial, meaning then number of elements in the zero_array

    zero_array: array containing the zeros of the polynomial

    Output: An array of N coefficients, with entry [0] being the coefficient of x^{n-1}
     of the polynomial with the zeros in zero_array.
    """
    coeffs = zeros(N, dtype=np.complex)
    for m in range(1, N+1):
        index_tracker = array([m - x for x in range(m)])
        coeff = 0
        first_iteration = True
        while True:
            if not first_iteration:
                index_tracker = get_next_index_tracker(N, index_tracker)
            else:
                first_iteration = False

            term = 1
            for i in range(m):
                term *= zero_array[index_tracker[i]-1]
            
            coeff += term

            if index_tracker[-1] == N - m + 1:
                break

        if m % 2 == 1:  # (-1)^n term
            coeff = -coeff

        coeffs[m-1] = coeff

    return coeffs


def f(t): # function f(t) defined in the paper. For my purposes I use f(t) = t
    return t


def f_p(t): # the derivative of f(t)
    return 1


def g(t, T):  # the function g(t) as defined in the paper.
    return f_p(t)/(f(T) - f(0))


def get_yp(n, degree, y, cms, yms, t, T=1): # returns a y'_n (0)
    """
    This create a single entry n of the system of differential equations vector F(t)
    It is a helper function for create_yp_vector.

    Input:

    n: The index of the entry of the vector we are considering

    y: The numerical approximation of vector F(t)

    cms: The array of coefficients c explained in the paper

    yms: The array of gamma values explained in the paper

    t: The time at which the vecotr F(t) is being approximated

    Output: entry n of F(t) at time t.
    """
    sum_term = 0+0j
    division_term = 1+0j
    for m in range(degree):
        if m != n:
            division_term *= (y[n] - y[m])
        sum_term += (cms[m] - yms[m]) * pow(y[n], degree - (m+1))
    return -g(t, T) * (sum_term / division_term)


def create_yp_vector(degree, y, cms, yms, t):
    """
    Creates the  F(t) vector at time t. This will be passed to our integrator.

    Input:

    degree: integer denoting the degree of the polynomial

    cms: The array of coefficients c explained in the paper

    yms: The array of gamma values explained in the paper

    t: The time at which the vector F(t) is being approximated

    Output: The vector F(t) at time t.
    """

    yp_vector = zeros(degree, dtype=complex)
    for n in range(degree):
        yp_vector[n] = get_yp(n, degree, y, cms, yms, t)
    return yp_vector



def runge_kutta4(f, x, t, h):
    """
    Runge Kutta integrator of fourth order. Same as used in the homework.
    """
    s1 = f(x, t)
    s2 = f(x + h/2 * s1, t + h/2)
    s3 = f(x + h/2 * s2, t + h/2)
    s4 = f(x + h * s3, t + h)
    return x + h/6 * (s1 + 2*s2 + 2*s3 + s4)


class Polynomial:  # This is a class that allows to store essential variables that are needed in the algorithm
    def __init__(self, degree, coeffs, bounds=[-100, 100], create_using_zeros=True):
        """
        degree:  The degree of the polynomial

        coeffs:  This can be one of two things. If create_using_zeros = True, then coeffs should be the zeros
        of the polynomial. If create_using zeros = False, then coeffs should be coefficients of the polynomial.

        bounds: array of 2 floats. These bounds that will be used to create an initial guess for the zeros. First
        entry is for the lower bound and second entry is for upper bound

        create_using_zeros: Specifies the way that the polynomial object is created. True implies it takes zeros
        and generates its coefficients. False implies it is given coefficients.
        """

        self.degree = degree
        self.coeffs = array(coeffs)
        if create_using_zeros: # Here if we specify the zeros of the polynomial, we generate the coefficients.
            self.coeffs = generate_coeffs(degree, self.coeffs)

        # We generate initial guesses based on specified bounds
        self.time_zeros = generate_initial_guesses(degree, bounds[0], bounds[1])

        # These variables will be used  to keep track of the time progression and zero progression.
        # Mainly used for being able to plot the results later.
        self.time_zero_progression = None # This will be an array that stores the zero progression over time
        self.time_progression = None # This will be an array that stores the time progression.

        # Debugging prints
        # print("Polynomial created with:")
        # print("coeffs = ", self.coeffs)
        # print("guesses = ", self.time_zeros)

    def print_properties(self):
        """ A simple debugging tool that prints some of the variable value stored in the object """
        print("degree:", self.degree)
        print("coeffs:", self.coeffs)
        print("time_zeros", self.time_zeros)
        print("yms", self.yms)

    def create_yp_vector(self, y, t):
        """ A wrapper which has only two parameters y,t allowing us to pass it to the  """
        return create_yp_vector(self.degree, y, self.coeffs, self.yms, t)

    def evaluate_time_zeros(self, h, T):
        """ This is a function that iterates the runge_kutta integrator.
        Input:

        h: timestep.

        T: the value until which our integration must go on. Explained in the paper.

        Output: None, but the time progression, time dependent zero progression and the final evaluated zeros are stored
        in the polynomial.
        """

        self.yms = generate_coeffs(self.degree, self.time_zeros) # We generate the gamma values
        times = [0] # We start at time = 0
        ys = [self.time_zeros]  # our F(0) is stored into the progression array ys
        while times[-1] < T :  # we iterate. This is a copy of the iterator used in homework 6.
            y_new = runge_kutta4(self.create_yp_vector, ys[-1], times[-1], h)
            ys.append(y_new)
            times.append(times[-1] + h)
        # We store the values in the polynomial.
        self.time_zeros = ys[-1]
        self.time_zero_progression = array(ys)
        self.time_progression = array(times)


def plot_zero_progression(p, correct_zeros):
    """
     Function used for plotting purposes. Take a polynomial and its correct zeros and plot them.
     Plots 'x' for where the progression starts, '^' for the correct zero.

    Input:

     p: polynomial whose time dependent zero progression we want to plot

    correct_zeros: The correct zeros of the polynomial

    Output: None. The zero progression is plotted. The correct zeros are plotted to see if the
    progression of the time zeros is heading in the right direction.
    """
    N = p.degree
    gca().set_color_cycle(None) # Used to reset colors, so that we can have the '^', the 'x' and line of the same color.
    for zero_index in range(N):
        plot(p.time_zero_progression.real[:, zero_index], p.time_zero_progression.imag[:, zero_index])
    gca().set_color_cycle(None)
    for zero_index in range(N):
        plot(correct_zeros.real[zero_index], correct_zeros.imag[zero_index], "^")
    gca().set_color_cycle(None)
    for zero_index in range(N):
        plot(p.time_zero_progression.real[0, zero_index], p.time_zero_progression.imag[0, zero_index], "x")
        print("Initial point {} let to the zero {}".format
              (p.time_zero_progression[0][zero_index], p.time_zero_progression[-1][zero_index]))
        print()

def evaluate_polynomial(N, coeffs, value):
    """
    Evaluates a polynomial at the given value

    Input:

    N: degree of polynomial

    coeffs: coefficient array, with coeffs[0] being the constant in front of x^{N-1} and coeffs[-1]
    being the constant term.

    value: The value at which we want to evaluate the polynomial

    Output: The value of the polynomial when z = value
    """
    s = 1 + 0j
    for i in range(N):
        s *= value
        s += coeffs[i]
    return s


if __name__ == "__main__": # Here I played around with the algorithm to see what results I could get.

    # Example where we expand from a set of zeros
    T= 1
    N = 4  # degree
    iteration_number = 10 #number of iterations of polynomial

    if __name__ == '__main__':
        zero_list = [6,-10, 3+5j, 6-4j]  # the zeros we want our polynomial to have

        # create the polynomial using these zeros
        # our initial guess for y_n(0) we say is bounded between -10 and 10 on both the
        # imaginary parts.
        p = Polynomial(N, zero_list, [-10, 10])

        for i in range(iteration_number):
              p.evaluate_time_zeros(0.01, T)  # Perform iterative algorithm

        print(p.time_zeros)  # print out the zeros we found

        # To create with coefficients instead of zeros, do the following:

        coeffs = [4, 5, 10, -3]

        p = Polynomial(N, coeffs, [-10, 10], False) # False to specify that we are generating polynomial from coeffs.
        print(p.coeffs)






