# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import interpolation as ip
import integration as ig

#----------------------------------#
#   TEST CUBIC SPLINES METHOD      #
#----------------------------------#

n = 100
f = lambda x: np.cos(x)
x = np.linspace(0, np.pi, n);
y = f(x);
y2 = ip.second_derivate(x, y, 0, 0, 999999);
mp.plot(x, y2);
mp.title("Second derivate of cosinus");
mp.show();

y_inter = np.zeros(n);
for i in range(0, n - 1):
    y_inter[i] = ip.cubic_spline_meth(x, y, y2, x[i]);

mp.plot(x, y_inter);
mp.title("Interpolation of cosinus with cubic_spline_meth");
mp.show();

#----------------------------------#
#     TEST INTEGRATION METHODS     #
#----------------------------------#

print "Rectangle integration method on cosinus between 0 and pi:"

print "Simpson method on cosinus between 0 and pi:"

print "Monte-Carlo method on cosinus between 0 and pi:"

print "Length of cosinus between 0 and pi:"
