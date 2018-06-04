''' Parameters used in test and benchmark methods '''
from __future__ import division, print_function, absolute_import

from random import random

from scipy.optimize import zeros as cc


def f1(x):
    return x*(x-1.)


def f2(x):
    return x**2 - 1


def f3(x):
    return x*(x-1.)*(x-2.)*(x-3.)


def f4(x):
    if x > 1:
        return 1.0 + .1*x
    if x < 1:
        return -1.0 + .1*x
    return 0


def f5(x):
    if x != 1:
        return 1.0/(1. - x)
    return 0


def f6(x):
    if x > 1:
        return random()
    elif x < 1:
        return -random()
    else:
        return 0


description = """
f2 is a symmetric parabola, x**2 - 1
f3 is a quartic polynomial with large hump in interval
f4 is step function with a discontinuity at 1
f5 is a hyperbola with vertical asymptote at 1
f6 has random values positive to left of 1, negative to right

of course these are not real problems. They just test how the
'good' solvers behave in bad circumstances where bisection is
really the best. A good solver should not be much worse than
bisection in such circumstance, while being faster for smooth
monotone sorts of functions.
"""

methods = [cc.bisect,cc.ridder,cc.brenth,cc.brentq]
mstrings = ['cc.bisect','cc.ridder','cc.brenth','cc.brentq']
functions = [f2,f3,f4,f5,f6]
fstrings = ['f2','f3','f4','f5','f6']
