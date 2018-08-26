from __future__ import division, print_function
from pylab import *


# Put here the functions ftcs, wave_first_step, wave
# (see below to figure out how they are used)
def ftcs(u, dx, dt, D):
    sigma = D/(dx*dx) * dt
    u_old = copy(u)
    for i in range(1, len(u)-1):
        u[i] = (1-2*sigma)*u_old[i] + sigma*(u_old[i+1] + u_old[i-1])
        print(u[i], u_old[i])
    return u

def wave_first_step():
    return

def wave():
    return


#--- Definition of the spatial grid -------------

Nmesh = 100
dx    = 1./Nmesh
x     = arange(0, 1.+dx/2., dx)

#--- Problem 1 ----------------------------------
D   = 1.e-3
#setting dt in such a way that sigma = 1/4
dt  = dx**2/(4.*D)
u = x*(1-x)*tanh(50*(0.5-x))
# u = zeros_like(x)
# u[Nmesh//4: 3*Nmesh//4] = 1.


Endtime = 3.
Nloops  = int(Endtime / dt)+1
hold(False)
for i in xrange(Nloops):
    u = ftcs(u, dx, dt, D)
    plot(x, u)
    title("time={:4.2f}".format(i*dt))
    draw()
    pause(0.001)


# #--- Problem 2 ----------------------------------
# c = 1.
# #setting dt in such a way that cfl = 1/4
# dt  = dx/(4*c)
#
# f = zeros_like(x)
# f[45:56] = sin(((x[45:56]-x[45])/(x[55]-x[45]) )*pi )**2
# g = zeros_like(x)
# g[1:-1] = c*(f[2:]-f[:-2])/(2*dx)
#
# uold, u = wave_first_step(f, g, dx, dt, c)
#
# Endtime = 1.
# Nloops  = int(Endtime / dt)+1
# hold(False)
# for i in xrange(1, Nloops):
#     plot(x, u)
#     ylim(-1, 1)
#     grid(True)
#     title("time={:4.2f}".format(i*dt))
#     draw()
#     pause(0.001)
#     uold, u = wave(uold, u, dx, dt, c)
#
#
#
# from IPython import embed
# embed()
