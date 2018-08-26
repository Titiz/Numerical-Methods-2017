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
    if f(a)*f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))

    iterations = 0
    while True:
        iterations +=1
        #print("a={},   b={},   b-a={}, b+a/2={}".format(repr(a), repr(b), repr(b-a), repr((b+a)/2)))
        midpoint = (a+b)/2.
        if fabs(f(midpoint)) < 1e-10:
            break
        sfm = sign(f(midpoint))
        #print("f(midpoint) = " +  repr(f(midpoint)))
        if sfm == sign(f(a)):
            a = midpoint
        elif sfm == sign(f(b)):
            b = midpoint
        elif sfm == 0:
            #lucky case: we found an exact zero!
            break
        else:
           raise ValueError("Something went horribly bad: sign(f(midpoint))={}\n".format(sfm))
    print(iterations)
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
    iterations = 0
    while True:
        iterations += 1
        time.sleep(0.5)
        newa = a - f(a)/fp(a)
        print(iterations)
        print("a={}, newa={},  |newa-a|={}".format(repr(a), repr(newa), repr(fabs(a-newa))))
        print("cos(a) =", repr(fp(a)))
        print("sin(a) =", repr(f(a)))
        print("sin(a)/cos(a) =", repr(f(a)/fp(a)))
        print("a - sin(a)/cos(a) =", newa)
        print("fabs(a-newa)/fabs(newa)= ", fabs(a-newa)/fabs(newa))
        if (isnan(fabs(a-newa)/fabs(newa))):
            break
        if (fabs(a-newa)/fabs(newa) < tol):
            break
        else:
            a = newa
    print(iterations)
    return newa



def regula_falsi(f, a, b, tol = 1e-10):
    if a > b:
        c=b
        b=a
        a=c 
    if a==b:
        raise ValueError("regula_falsi called with a==b\n")
    if f(a)*f(b) >= 0.:
        raise ValueError("The interval does not bracket a zero! f({})={}; f({})={}\n".format(a, f(a), b, f(b)))
    
    iterations = 0
    while True:
        iterations += 1
        c = (b*f(a) - a*f(b)) / (f(a) - f(b))
        if (f(c) == 0):
            break
        if f(a)*f(c) < 0:
            b = c
        else:
            a = c
    print(iterations)
    return c

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
        return sin(x-1)

    def fp(x):
        return cos(x-1)

    #new2 = Newton(lambda x: x*x*x*x + x*x*x, lambda x: 4*x*x*x + 3*x*x, 0.5)

    #Problem 4:

    #print(regula_falsi(f, -1, 0.5))
    #print(bisection(f, -1, 0.5))
    #print(Newton(f, fp, -1))


    #Problem 5:


    print(Newton(f, fp, 4.))
    print(Newton(f, fp, 2.))

