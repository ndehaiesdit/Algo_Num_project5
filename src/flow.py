# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import interpolation as ip
import integration as ig
import charg_f as cf

# Tracement de l'aile

xup, yup, xdown, ydown = cf.load_foil("boe103.dat")
y2up = ip.second_derivate(xup, yup, 0, 0, 999999)
hup = np.max(yup)
y2down = ip.second_derivate(xdown, ydown, 0, 0, 999999)
hdown = np.min(ydown)

# Retourne la fonction qui associe l'ordonnée de l'aile haute à l'abscisse x
def y_wing_up(x):
    return lambda x: ip.cubic_spline_meth(xup, yup, y2up, x)

# Retourne la fonction qui associe l'ordonnée de l'aile basse à l'abscisse x
def y_wing_down(x):
    return lambda x: ip.cubic_spline_meth(xdown, ydown, y2down, x)

# Tracé des lignes de courants

# Retourne les ordonnées de la ligne de champ connaissant les caractéristiques h et y_wing de l'aile
# et la valeur lamb de la ligne de champ.
def airflow(x, y_wing, lamb, h):
    return (1. - lamb) * y_wing(x) + lamb*3.*h

def derivate(f, epsilon):
    return lambda x: (f(x + epsilon) - f(x)) / epsilon

# Renvoie lambda en donnant une abscisse x, une ordonnée y, l'ordonnée du point de l'aile y_wing,
# et l'extremum h du point de l'aile
def reciprocal(x, y, y_wing, h):
    return (y - y_wing(x)) / (3 * h - y_wing(x))

# Renvoie la longueur de la ligne de champ rencontrée au point (x, y) entre le 
# connaissant les paramètres de l'aile y_wing et h
def length_airflow(x, y):
    h = hup
    lamb = reciprocal(x, y, y_wing_up, h);
    y_wing = y_wing_up
    if lamb < 0 or lamb > 1:
        h = hdown
        lamb = reciprocal(x, y, y_wing_down, h);
        y_wing = y_wing_down
    if lamb < 0 or lamb > 1:
        lamb = 1.
    
    f = lambda x: airflow(x, y_wing, lamb, h)
    return length(ig.simpson_meth, 100, derivate(f, 0.0001), xup(0), xup(xup.size - 1))

print length_airflow(0.2, 0.2) 

         
