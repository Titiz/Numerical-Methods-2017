from __future__ import division, print_function
from pylab import *


def logistic(x, t=None):
    """Logistic growth model.
Inputs
x:  the concentration of the population expressed as individuals/area
t:  dummy variable that allows to use this function with integrators written
    for non-autonomous problems

Returns the growth rate of the population (it can be negative)
"""
    # I know it's ugly to bolt-in the constants in this way.  However,
    # more general, elegant and useful methods for specifying a
    # function's parameters, based on object-oriented programming or on
    # closures may confuse those of you who are inexperienced
    # programmers.  If you know what I'm saying, feel free to use those
    # techniques. If you have no clue of what I'm talking about, and
    # would like to know, feel free to ask!
    k = 1.
    r = 1.
    x = float(x)
    return r * x * (1 - x / k)


def ho(x, t=None):
    """Harmonic oscillator.
Inputs
x:  a two-elements vector containing position and velocity (in this
    order) of the mass
t:  dummy variable that allows to use this function with integrators written
    for non-autonomous problems

Returns the right-hand side of equations (6.42) in the book
"""
    k = 1.  # ratio of spring constant and mass of the object
    t = float(t)
    x = array(x)
    if x.shape != (2,):
        raise ValueError("The state of the oscillator is described by a vector with two elements")
    return array([x[1],
                  - x[0]])


def dfp(x, t):
    """Damped forced pendulum.
Inputs
x:  a two-elements vector containing angular position and velocity (in this
    order) of the pendulum
t:  time

Returns the right-hand side of equations (6.42) in the book
"""
    k = 1.  # ratio of gravity acceleration to length of pendulum
    d = 0.1  # damping coefficient
    A = 1.5  # amplitude of forcing
    t = float(t)
    x = array(x)
    if x.shape != (2,):
        raise ValueError("The state of the pendulum is described by a vector with two elements")
    return array([x[1],
                  - k * sin(x[0]) - d * x[1] + A * sin(t)])


def Lorenz63(x, t=None):
    """Saltzman's equations studied by Lorenz in his 1963 paper on chaos.
Inputs
x:  a 3-elements vector
t:  dummy variable that allows to use this function with integrators written
    for non-autonomous problems

Returns the right-hand side of equations (6.53) in the book, but with
different parameters.

    """
    s = 10.  # Prandtl's number
    r = 20.  # 13.92655741 #16.    #28.   #Rayleigh number
    b = 8. / 3.  # aspect ratio
    x = array(x)
    if x.shape != (3,):
        raise ValueError("The state of Lorenz's system is described by a vector with 3 elements")
    return array([- s * (x[0] - x[1]),
                  - x[0] * (x[2] - r) - x[1],
                  x[0] * x[1] - b * x[2]])


# --------------------------------------------------------------------------

def Euler(f, x, t, h):
    """Performs one step of the Euler integration.

Inputs
x:  the state vector at time t
t:  the current time
h:  the time step
f:  a function that returns a vector with the same dimensions of x when
    called as  f(x,t)

Returns the estimate of x at time t+h computed with Euler's method.
"""
    return x + h * f(x, t)


def trapezoid(f, x, t, h):
    return x + h / 2 * (f(x, t) + f(x + h * f(x, t), t + h))


def runge_kutta4(f, x, t, h):
    s1 = f(x, t)
    s2 = f(x + h / 2 * s1, t + h / 2)
    s3 = f(x + h / 2 * s2, t + h / 2)
    s4 = f(x + h * s3, t + h)
    return x + h / 6 * (s1 + 2 * s2 + 2 * s3 + s4)


def iterate_solver(sol, fun, xstart, tstart, h, tend):
    """Iterate a one-step solver with time step h starting from initial
condition xstart at tstart, until tend is reached.

Inputs
sol:    a function that implements a one-step solver. sol(f,x,t,h) must
        yield an estimate of x(t+h)
fun:    the vector field that defines the (system) of ordinary differential
        equation(s). The time derivative of the state x must be given by f(x,t)
xstart: the initial condition, a vector of n elements, or a float for 1D
        problems.
tstart: the time corresponding to the initial condition. It's ignored for
        autonomous problems (that is, problems where f does not explicitly
        depend on time).
tend:   the time after which the iterations must stop.

Returns
xs:    an (m+1)xn matrix containing the state vectors of the solution.
       xs[0,:] is equal to the initial condition xstart. m is the number of
       time steps performed. n is the number of elements of xstart.
times: the times correspoding to the vectors in xs. times[0] is tstart
"""
    xs = [xstart]
    times = [tstart]
    iter = 0
    while times[-1] < tend:
        iter += 1
        xnew = sol(fun, xs[-1], times[-1], h)
        xs.append(xnew)
        times.append(tstart + iter * h)

    return array(xs), array(times)


# ----------------------------------------------


# --- Exercise 2


if __name__ == "__main__":

        figure()

        x, t = iterate_solver(Euler, logistic, 0.01, 0, 0.1, 15.)
        plot(t, x)

        x, t = iterate_solver(trapezoid, logistic, 0.01, 0, 0.1, 15.)
        plot(t, x)

        x, t = iterate_solver(runge_kutta4, logistic, 0.01, 0, 0.1, 15.)
        plot(t, x)

        legend(["Euler", "Trapezoid", "RK4"], loc="upper left")

        figure()

        x, t = iterate_solver(Euler, ho, [1., 0.], 0, 0.05, 30.)
        plot(t, x[:, 0])

        x, t = iterate_solver(trapezoid, ho, [1., 0.], 0, 0.05, 30.)
        plot(t, x[:, 0])

        x, t = iterate_solver(runge_kutta4, ho, [1., 0.], 0, 0.05, 30.)
        plot(t, x[:, 0])

        legend(["Euler", "Trapezoid", "RK4"], loc="upper left")

        figure()

        x, t = iterate_solver(Euler, dfp, [1., 0.], 0, 0.01, 100.)
        plot(t, x[:, 0])

        x, t = iterate_solver(trapezoid, dfp, [1., 0.], 0, 0.01, 100.)
        plot(t, x[:, 0])

        x, t = iterate_solver(runge_kutta4, dfp, [1., 0.], 0, 0.01, 100.)
        plot(t, x[:, 0])
        legend(["Euler", "Trapezoid", "RK4"], loc="upper left")
        plt.show()



