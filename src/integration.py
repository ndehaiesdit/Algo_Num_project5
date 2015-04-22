# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr

#----------------------------------#
#         INTEGRAL METHODS         #
#----------------------------------#

def rect_meth(x, y, n, f):
    a = (y-x)/n
    s = 0.0
    t = np.arange(0, n)
    for i in range(0,n):
        s += f(x + i*a)
    return a*s
    
def simpson_meth(x, y, n, f):
    a = (y-x)/n
    s1 = 0.0
    s2 = 0.0
    t = np.arange(1, n/2)
    u = np.arange(1, (n/2)+1)
    for i in t:
        s1 += f(x + 2*i*a)
    s1 = 2*s1
    for i in u:
        s2 += f(x + (2*i-1)*a)
    s2 = 4*s2
    return (a/3)*(f(x) + s1 + s2 + f(y))

def monte_carlo_meth(x, y, n, f):
    max = 0.0
    s = 0.0
    for i in range(n+1):
        a = f(x + y*float(i+1)/n)
        if (max < a):
            max = a
    alpha = (b-a)*max
    for i in range(n):
        x_random = (b-a)*npr.rand() + x
        y_random = max*npr.rand()
        if (y_random <= max):
            s += (alpha/n) * f(x_random)
    return s


#--------------------#
#       LENGTH       #
#--------------------#

def length(I, n, df, x, y):
    length = lambda x: np.sqrt(1 + (df(x)*df(x)))
    return I(x, y, n, length)

def coef_to_fonction(A, B, C, D):
    return length x: A*(x**3) + B*(x**2) + C*x + D

def coef_to_deriv_fonction(A, B, C, D):
    return length x: A*3*(x**2) + B*2*x + C    


