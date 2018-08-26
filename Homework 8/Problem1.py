from __future__ import division, print_function
from pylab import *


# Put here the functions ftcs, wave_first_step, wave
# (see below to figure out how they are used)
def ftcs(u, dx, dt, D):
    sigma = D/(dx*dx) * dt # Defines the sigma value
    u_next = zeros_like(u) # Create an empty vector
    for i in range(1, len(u)-1): # iterate and create the u vector at the next t-value.
        u_next[i] = (1-2*sigma)*u[i] + sigma*(u[i+1] + u[i-1])
    return u_next  # return the u vector at the next x-coord.

#--- Definition of the spatial grid -------------
Nmesh = 100
dx    = 1./Nmesh
x     = arange(0, 1.+dx/2., dx)

if __name__ == "__main__":
    #--- Problem 1 ----------------------------------
    D   = 1.e-3
    #setting dt in such a way that sigma = 1/4
    dt  = dx**2/(4.*D)
    u = x*(1-x)*tanh(50*(0.5-x))
    # u = zeros_like(x)
    # u[Nmesh//4: 3*Nmesh//4] = 1.


    Endtime = 500.
    Nloops  = int(Endtime / dt)+1
    hold(False)
    for i in xrange(Nloops):
        u = ftcs(u, dx, dt, D)
        plot(x, u)
        title("time={:4.2f}".format(i*dt))
        draw()
        pause(0.001)
    show()


