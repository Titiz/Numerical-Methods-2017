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


if __name__ == "__main__":
    # --- Exercise 4
    z_min = 18.0
    z_max = 19.0
    x_min = 0.2
    x_max = 1.2
    Ndots = 201
    x_space = linspace(x_min, x_max, Ndots)
    z_space = linspace(z_min, z_max, Ndots)
    results = zeros((Ndots, Ndots))
    type1s = [] #  used for file output debugging
    type2s = [] #  used for file output debugging
    for x_index in range(Ndots):
        for z_index in range(Ndots):

            print("testing with values x = {},  z = {}".format(x_space[x_index], z_space[z_index]))

            x1, t1 = iterate_solver(runge_kutta4, Lorenz63,
                                    [x_space[x_index], -x_space[x_index], z_space[z_index]],
                                    0, 0.01, 60)

            last_x = x1[-1, 0]
            last_y = x1[-1, 1]
            last_z = x1[-1, 2]
            print("x at time", t1[-1], "is", last_x)
            print("y at time", t1[-1], "is", last_y)
            print("z at time", t1[-1], "is", last_z)
            if last_x < 0:  # Type 1 converges to (-7.118, -7.118, 19)
                results[x_index, z_index] = -1
                print("starting condition is of type 1 \n")
                type1s.append([last_x, last_y, last_z])  #  used for file output debugging
            elif last_x > 0:  # Type 2 converges to (7.118, 7.118, 19)
                results[x_index, z_index] = 1
                print("starting condition is of type 2 \n")
                type2s.append([last_x, last_y, last_z])  #  used for file output debugging
            #   # Graphs used for manual observation of the end conditions debugging purposes
            # figure()
            # plot(t1, x1[:, 0])
            # plot(t1, x1[:, 1])
            # plot(t1, x1[:, 2])
            # plot(x1[:, 2], x1[:, 0])
            # plot(x1[:, 2], x1[:, 1])
            # plot(x1[:, 0], x1[:, 1])
            # show()

    X, Z = meshgrid(x_space, z_space)
    color_mesh = pcolormesh(X, Z, results)
    axis('image')
    plt.yticks(np.arange(18, 19, 1 / (Ndots-1)))
    plt.xticks(np.arange(0.2, 1.2, 1 / (Ndots - 1)))
    plt.xlabel("x value")
    plt.ylabel('z value')
    grid()
    plot()
    show()

    ''' Extra attempt:
    Any three points in the form (0, 0, z) for any real z will converge to (0,0,0) as t goes to infinity '''

    #  # Exploration of the two types
    #  # Type 1:
    # x1, t1 = iterate_solver(runge_kutta4, Lorenz63,
    #                         [-1, 1, 18.5],
    #                         0, 0.01, 100)
    # figure()
    # plot(t1, x1[:, 0])  #  x vs time
    # plot(t1, x1[:, 1])    #  y vs time
    # plot(t1, x1[:, 2])    #  z vs time
    # plot(x1[:, 2], x1[:, 0])   #  z vs x
    # plot(x1[:, 2], x1[:, 1])   #  z vs y
    # plot(x1[:, 0], x1[:, 1])   #  x vs y
    # print(x1[-1, 0])
    # print(x1[-1, 1])
    # print(x1[-1, 2])
    #   # Type 2:

    # x1, t1 = iterate_solver(runge_kutta4, Lorenz63,
    #                         [-1, 1, 18.4],
    #                         0, 0.01, 100)
    # figure()
    #
    # plot(t1, x1[:, 0])
    # plot(t1, x1[:, 1])
    # plot(t1, x1[:, 2])
    # plot(x1[:, 2], x1[:, 0])
    # plot(x1[:, 2], x1[:, 1])
    # plot(x1[:, 0], x1[:, 1])
    # print(x1[-1, 0])
    # print(x1[-1, 1])
    # print(x1[-1, 2])
    # show()

    #  Output used for tracking all of the x,y,z coordinates at the last time for type 1 and type 2 outputs.
    f = open("type1s.txt", "w")
    for element in type1s:
        f.write("{}, {}, {}".format(element[0], element[1], element[2]) + "\n")
    f.close()
    f = open("type2s.txt", "w")
    for element in type2s:
        f.write("{}, {}, {}".format(element[0], element[1], element[2]) + "\n")
    f.close()

# ---------------------------------

