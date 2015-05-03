# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import interpolation as ip
import integration as ig
import charg_f as cf

# Tracement de l'aile

xhaut, yhaut, xbas, ybas = cf.load_foil("boe103.dat")

y2haut = ip.second_derivate(xhaut, yhaut, 0, 0, 999999)
hhaut = np.max(yhaut)

n = np.size(xhaut)

y_interhaut = np.zeros(n)

for i in range(0, n - 1):
    y_interhaut[i] = ip.cubic_spline_meth(xhaut, yhaut, y2haut, xhaut[i])

y2bas = ip.second_derivate(xbas, ybas, 0, 0, 999999)
hbas = np.min(ybas)

n = np.size(xbas)

y_interbas = np.zeros(n)

for i in range(0, n - 1):
    y_interbas[i] = ip.cubic_spline_meth(xbas, ybas, y2bas, xbas[i])

mp.plot(xhaut, y_interhaut, 'g')
mp.plot(xbas, y_interbas, 'r')
mp.gca().set_aspect('equal', adjustable='box')
mp.show()    


#Tracé des lignes de courants

def ligne(x, y, lamb, h):
    n = np.size(x)
    y_ligne = np.zeros(n)

    for i in range(0, n):
        y_ligne[i] = (1 - lamb) * y[i] + lamb*3*h

    return y_ligne

def plot_isobar(x, y, y_extremum, dy):
    lamb = 0.
    dlamb = np.abs(dy / (3. * y_extremum))
    while lamb < 1:
        mp.plot(x, ligne(x, y, lamb, y_extremum), 'black')
        lamb = lamb + dlamb



def derivate(f, epsilon):
    return lambda x: (f(x + epsilon) - f(x)) / epsilon

# Renvoie lambda en donnant une ordonnée y, l'ordonnée du point de l'aile y_wing,
# et l'extremum h du point de l'aile
def reciprocal(y, y_wing, h):
    return (y - y_wing) / (3 * h - y_wing)



#def pressure(x, y):
#    line = lambda a: ligne(xhaut, y, 
#    return ig.length(ig.mont_carlo_meth, 99, derivate(line, 0.01), 0, x) - x

plot_isobar(xhaut, y_interhaut, hhaut, 0.01)
plot_isobar(xbas, y_interbas, hbas, 0.01)
mp.title('Airflow arround the airfoil of a Boeing')
mp.gca().set_aspect('equal', adjustable='box')
mp.show() 

         
