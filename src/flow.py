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

# Retourne l'ordonnée de l'aile haute à l'abscisse x
def y_wing_up(x):
    return ip.cubic_spline_meth(xup, yup, y2up, x)

# Retourne l'ordonnée de l'aile basse à l'abscisse x
def y_wing_down(x):
    return ip.cubic_spline_meth(xdown, ydown, y2down, x)

# Retourne les ordonnées de la ligne de champ 
# connaissant les caractéristiques h et y_wing de l'aile
# et la valeur lamb de la ligne de champ.
def airflow(x, y_wing, lamb, h):
    return (1. - lamb) * y_wing(x) + 3. * lamb * h

# Retourne une approximation de la dérivée de f à partir d'un coefficient epsilon
def derivate(f, epsilon):
    return lambda x: (f(x + epsilon) - f(x)) / epsilon

# Renvoie lambda en donnant une abscisse x, une ordonnée y, 
# l'ordonnée du point de l'aile y_wing, et l'extremum h du point de l'aile
def reciprocal(x, y, y_wing, h):
    y_wing_x = y_wing(x)
    return (y - y_wing_x) / (3 * h - y_wing_x)

# Renvoie la longueur de la ligne de champ rencontrée au point (x, y) entre le 
# connaissant les paramètres de l'aile y_wing et h
def length_airflow(x, y, iterations):
    h = hup
    lamb = reciprocal(x, y, y_wing_up, h);
    y_wing = y_wing_up
    if lamb < 0 or lamb > 1:
        h = hdown
        lamb = reciprocal(x, y, y_wing_down, h);
        y_wing = y_wing_down
    if lamb < 0 or lamb > 1:
        return xup[xup.size - 1] - xup[0]
    
    f = lambda x: airflow(x, y_wing, lamb, h)
    return ig.length(ig.simpson_meth, iterations, derivate(f, 0.0001), 
                     xup[0], xup[xup.size - 1])

# Plot les lignes de courant du coté de l'aile y_wing avec h son ordonnée extremum,
# dy le pas séparant les lignes et n le nombre de point des lignes
def plot_airflows_side(dy, n, y_wing, h):
    x_airflow = np.linspace(xup[0], xup[xup.size - 1], n)
    y_airflow = np.zeros(n)
    lamb = 0.
    dlamb = np.abs(dy / (3. * h))
    while lamb < 1:
        for i in range (0, n):
            y_airflow[i] = airflow(x_airflow[i], y_wing, lamb, h)
        mp.plot(x_airflow, y_airflow, color='black')
        lamb += dlamb

# Affiche les lignes de courant :
# dy est le pas séparant les lignes et n le nombre de point des lignes
def display_airflows(dy, n):
    plot_airflows_side(dy, n, y_wing_up, hup)
    plot_airflows_side(dy, n, y_wing_down, hdown)
    mp.gca().set_aspect('equal', adjustable='box')
    mp.title("Airflows")
    mp.show();
    

# Retourne la matrice de contenant les longueurs des lignes 
# de champ avec matrix[0][0] = length_airflow(x, y) 
# et un pas horyzontal et vertical step
# iterations est le nombre d'itérations pour le calcul des longueurs
def matrix_length(x, y, width, height, step, iterations) :
    n = int(height / step)
    m = int(width / step)
    matrix = np.zeros((n, m))
    for i in range(0, n):
        for j in range(0, m):
            matrix[i,j] = length_airflow(x, y, iterations);
            x += step;
        x = 0;
        y -= step;
    return matrix

         
