from __future__ import division, print_function
from pylab import *


# Put here the functions ftcs, wave_first_step, wave
# (see below to figure out how they are used)

def wave_first_step(f, g, dx, dt, c):
    sigma_2 = c * c * dt * dt / (dx * dx)  # create the sigma^2 value
    u = zeros_like(f)  #create an empty vector u
    for i in range(1, len(f)-1): # fill vector u up
        u[i] = (1-sigma_2) * f[i] + \
            sigma_2/2 * (f[i+1] + f[i-1]) - \
            dt * g[i]
    return array([f, u])  # return f as oldu and u as u.


def wave(uold, u, dx, dt, c):
    sigma_2 = c * c * dt * dt / (dx * dx) # create the sigma^2 value
    u_next = zeros_like(u) # create an empty vector at next t-value
    for i in range(1, len(u)-1): # iterate and fill the vector at next t-value
        u_next[i] = 2*(1-sigma_2)*u[i] + sigma_2*(u[i+1] + u[i-1]) - uold[i]
    return array([u, u_next])  #return u as oldu and u_next as u.


#--- Definition of the spatial grid -------------
Nmesh = 100
dx    = 1./Nmesh
x     = arange(0, 1.+dx/2., dx)


#--- Problem 2 ----------------------------------

if __name__ == "__main__":
    c = 1.
    #setting dt in such a way that cfl = 1/4
    dt = dx/(4*c)
    f = zeros_like(x)
    f[45:56] = sin(((x[45:56]-x[45])/(x[55]-x[45]) )*pi )**2
    g = zeros_like(x)
    g[1:-1] = c*(f[2:]-f[:-2])/(2*dx)

    uold, u = wave_first_step(f, g, dx, dt, c)

    Endtime = 1.
    Nloops  = int(Endtime / dt)+1
    hold(False)
    for i in xrange(1, Nloops):
        plot(x, u)
        ylim(-1, 1)
        grid(True)
        title("time={:4.2f}".format(i*dt))
        draw()
        pause(0.001)
        uold, u = wave(uold, u, dx, dt, c)



