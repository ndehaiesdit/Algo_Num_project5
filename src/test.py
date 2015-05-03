# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import math
import numpy.random as npr
import interpolation as ip
import integration as ig
import charg_f as cf

#----------------------------------#
#   TEST CUBIC SPLINES METHOD      #
#----------------------------------#

n = 100
f = lambda x: np.cos(x)
x = np.linspace(0, np.pi, n);
y = f(x);
y2 = ip.second_derivate(x, y, 0, 0, 999999);
y2_real = -np.cos(x)
mp.plot(x, y2, label='Calculated second derivate', color='red');
mp.plot(x, y2_real, label='Real second derivate', color='blue');
mp.title("Second derivate of cosinus");
mp.legend(loc='upper left')
mp.show();

y_inter = np.zeros(n);
for i in range(0, n - 1):
    y_inter[i] = ip.cubic_spline_meth(x, y, y2, x[i]);

mp.plot(x, y_inter, label='Interpolation of cosinus', color='red');
mp.plot(x, y, label='True cosinus', color='blue')
mp.title("Interpolation of cosinus with cubic_spline_meth");
mp.legend(loc='upper right')
mp.show();

#----------------------------------#
#     TEST INTEGRATION METHODS     #
#----------------------------------#

print "Rectangle integration method on cosinus between 0 and pi/2 with 100 iterations:"
print ig.rect_meth(0, np.pi/2., 100, f)

print "Simpson method on cosinus between 0 and pi/2 with 100 iterations:"
print ig.simpson_meth(0.0, np.pi/2., 100, f)

print "Middle point method on cosinus between 0 and pi/2 with 100 iterations:"
print ig.midpoint_meth(0.0, np.pi/2., 100, f)

#----------------------------------#
#     TEST LENGTH COMPUTING        #
#----------------------------------#

dfi = lambda x: cf.general_function(x,cf.ix,ig.derivate(cf.Poli))
dfe = lambda x: cf.general_function(x,cf.ex,ig.derivate(cf.Pole))

print "Length of cosinus between 0 and pi/2 with 100 iterations:"
print ig.length(ig.rect_meth,100, lambda x: -np.sin(x), 0.0, np.pi/2.)
print "Intrados length:"
print ig.length(ig.rect_meth,100,lambda x: cf.general_function(x,cf.ix,ig.derivate(cf.Poli)),0.0,1.0)
print "Extrados length:"
print ig.length(ig.rect_meth,100,lambda x: cf.general_function(x,cf.ex,ig.derivate(cf.Pole)),0.0,1.0)

ig.comp_conv(lambda x: np.sqrt(1+(dfi(x)*dfi(x))),0.0,1.0, 1.01639654139)
ig.comp_conv(lambda x: np.sqrt(1+(dfe(x)*dfe(x))),0.0,1.0, 1.078732871)
ig.comp_conv(np.cos, 0.0, np.pi/2., 1.0)
