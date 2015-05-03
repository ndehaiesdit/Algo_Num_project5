# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import charg_f as cf

#----------------------------------#
#         INTEGRAL METHODS         #
#----------------------------------#

def rect_meth(x, y, n, f):
    a = (y-x)/n
    s = 0.0
    t = np.arange(0, n)
    for i in t:
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

def midpoint_meth(x,y,n,f):
    s = 0.0
    a = (y-x)/n
    for k in np.arange(0,n):
        s = s + f(x+(k*a)+(a/2.))
    return a*s

def comp_conv(f,a,b,exact_value):
    nmax = 100
    
    x = np.arange(1.0, nmax, 1.0)

    t = np.arange(1.0, nmax, 1.0)
    u = np.arange(1.0, nmax, 1.0)
    v = np.arange(1.0, nmax, 1.0)
    w = np.arange(1.0, nmax, 1.0)

    for i in np.arange(0.0,t.size):
        t[i] = rect_meth(a,b,x[i],f)

    for i in np.arange(0.0,u.size):
        u[i] = simpson_meth(a,b,x[i],f)

    for i in np.arange(0.0,v.size):
        v[i] = midpoint_meth(a,b,x[i],f)

    for i in np.arange(0.0,v.size):
        w[i] = exact_value

    mp.clf()
    mp.semilogx()
    mp.plot(x, t, linewidth=1.0, label='Rectangle')
    mp.plot(x, u, linewidth=1.0, label='Simpson')
    mp.plot(x, v, linewidth=1.0, label='Midpoint')
    mp.plot(x, w, linewidth=1.0, label='Value of the integral')
    
    
    mp.title("Illustration of the convergence of the three integration methods");
    mp.xlabel('Number of subdivision points')
    mp.ylabel('Integral of f between a and b')
    mp.legend(loc='upper right')
    mp.show()

#--------------------#
#       LENGTH       #
#--------------------#

def length(I, n, df, x, y):
    length = lambda x: np.sqrt(1 + (df(x)*df(x)))
    return I(x, y, n, length)

# def coef_to_fonction(A, B, C, D):
#     return lambda x: A*(x**3) + B*(x**2) + C*x + D

# def coef_to_deriv_fonction(A, B, C, D):
#     return lambda x: A*3*(x**2) + B*2*x + C   
    

def derivate(Pol):
    n = Pol.shape[0]
    m = Pol.shape[1]
    dPol = np.zeros([n,m])
    
    for i in np.arange(0,n):
        dPol[i][0] = 0
        dPol[i][1] = 3*Pol[i][0]
        dPol[i][2] = 2*Pol[i][1]
        dPol[i][3] = Pol[i][2]
        
    return dPol
