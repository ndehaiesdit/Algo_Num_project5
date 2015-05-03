# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr

#----------------------------------#
#        CUBIC SPLINES METHOD      #
#----------------------------------#

# Approximation of the second derivate of a function taking a value y on each x,
# a derivate yp1 on the first x, a derivate ypn on the last.
# max is supposed to be a "huge" value in order to prevent a divergence
def second_derivate(x, y, yp1, ypn, max):
    n = np.size(x) - 1;
    u = np.zeros(n);
    y2 = np.zeros(n + 1);
    
    if(yp1 > max):
        y2[0] = u[0] = 0.0;
    else :
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[2]) - yp1);

    for i in range(1, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;

    if(ypn > max):
        qn = un = 0.0;
    else:
        qn = 0.5;
        un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);

    i = n - 1;
    while(i >= 0):
        y2[i] = y2[i] * y2[i + 1] + u[i];
        i = i - 1;

    return y2;

# The cubic spline interpolation method of a function f 
# which takes the value y for each x, a second derivate y2 for each x. 
# It returns f(value)
def cubic_spline_meth(x, y, y2, value):
    n = np.size(x) - 1;
    klo = 0;
    khi = n;

    while(khi - klo > 1):
        k = (khi + klo) >> 1;
        if(x[k] > value):
            khi = k;
        else:
            klo = k;

    h = x[khi] - x[klo];
    a = (x[khi] - value) / h;
    b = (value - x[klo]) / h;
    
    return a * y[klo] + b * y[khi] + ((a**3 - a) * y2[klo] + (b**3 - b) * y2[khi]) * h**2 / 6.0;
