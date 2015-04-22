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
x = np.linspace(0, np.pi, n);
y = np.cos(x);
y2 = ip.second_derivate(x, y, 0, 0, 999999);
mp.plot(x, y2);
mp.title("Derivee seconde de cosinus avec second_derivate");
mp.show();

y_inter = np.zeros(n);
for i in range(0, n - 1):
    y_inter[i] = ip.cubic_spline_meth(x, y, y2, x[i]);

mp.plot(x, y_inter);
mp.title("Interpolation de cosinus avec cubic_spline_meth");
mp.show();

#----------------------------------#
#     TEST INTEGRATION METHODS     #
#----------------------------------#

