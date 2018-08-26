from __future__ import division, print_function
from pylab import *


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


def bisection(f, a, b, tol = 1e-10):
    """Finds a zero of f by the bisection method

The sign of f(a) and f(b) must be opposite, and f must be a continuous function.
Input:
   f   function with only one float argument, that returns a float
   a   float, one end of the bracketing interval
   b   float, the other end of the bracketing interval

 tol   stop the iteration when (b-a)/2 < tol

Output:
   x such that f(x) is approximately 0.
"""
    # sanity checks to cover us from user distraction (n.b. usually user=me)
    if a > b:
        c=b
        b=a
        a=c
    if a==b:
        raise ValueError("Bisection called with a==b\n")
    if f(a)*f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    while (b-a) > tol:
        midpoint = (a+b)/2.
        sfm = sign(f(midpoint))
        if sfm == sign(f(a)):
            a = midpoint
        elif sfm == sign(f(b)):
            b = midpoint
        elif sfm == 0:
            #lucky case: we found an exact zero!
            return midpoint
        else:
           raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfm))
    return (a+b)/2.


def system(x, t=None):
    x = array(x)
    if x.shape != (2,):
        raise ValueError("The state of the oscillator is described by a vector with two elements")
    return array([x[1],
                  sin(x[1])])
    # return array([x[1], 4*x[0]])


def function(x):
    y, t = iterate_solver(runge_kutta4, system, [1, x], 0, h, 1)
    return y[-1, 0] + 1


if __name__ == "__main__":
    # Problem 1 #
    h = 0.01
    starting_v = bisection(function, -2, 2)
    xlabel("$t$", fontsize=20)
    ylabel("$y$", fontsize=20)
    y, t = iterate_solver(runge_kutta4, system, [1, starting_v], 0, h, 1)
    plot(t, y[:, 0])
    show()
