#the following import makes division work as in Python 3. Thus 2/3=0.66666....
#rather than 2/3=0 (in Python 2 division between integers yields an integer
#Python 3 has the specialized division operator // for integer division)
from __future__ import division 
from pylab import *
import time

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
    #if f(a)*f(b) >= 0.:
        #raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    while (2*) > tol:
        print("a={},   b={},   b-a={}, b+a/2={}".format(repr(a), repr(b), repr(b-a), repr((b+a)/2)))
        midpoint = (a+b)/2.
        sfm = sign(f(midpoint))
        print("f(midpoint) = " +  repr(f(midpoint)))
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

       

def Newton(f, fp, a, tol=1e-10):
    """Finds a zero of f by Newton's method

Input:
   f   function with only one float argument, that returns a float
  fp   the first derivative of f
   a   float

 tol   stop the iteration when |newa-a| < tol, where a, newa are extimates of the zero at consecutive iterations.

Output:
   x such that f(x) is approximately 0.
"""
    while True:
        newa = a - f(a)/fp(a)
        print("a={}, newa={},  |newa-a|={}".format(repr(a), repr(newa), repr(fabs(a-newa))))
        print("fp(a) =", repr(fp(a)))
        print("f(a) =", repr(f(a)))
        if isnan(newa):
            print("Probably tried to approach an assymptote.")
            if (a) >= 0:
                newa = float("Inf")
            elif (a) < 0:
                newa = -float("Inf")
            break
        if fabs(a-newa) < tol:
            break
        else:
            a = newa
    return newa




#------------------------
if __name__ == "__main__":
    print("Program started")
    #Problem 1:
    #print(bisection(sin, 1.e6*pi, (1.e6+1)*pi))

    # Problem 2:
    #bs = bisection(lambda x: (x-pi)*(x-pi), 3.1, 3.2)
    #print(repr(bs), ":" ,repr((bs-pi)*(bs-pi)))
    #new = Newton(lambda x: (x-pi)*(x-pi), lambda x: 2*x -2*pi, 1000000)
    #print(repr(new), ":", repr(((new-pi)*(new-pi))))

    #Problem 3:
    def f(x):
        return x*exp(-x)

    def fp(x):
        return exp(-x) + x * (-exp(-x))

    def g(x):
        return exp(-x)

    def gp(x):
        return -exp(x)
    
    new = Newton(f, fp, 2)
    print(new)
    #new2 = Newton(lambda x: x*x*x*x + x*x*x, lambda x: 4*x*x*x + 3*x*x, 0.5)