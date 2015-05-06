import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import interpolation as ip
import integration as ig
import charg_f as cf
import flow.py

def matrix_tranfo(n, p) :
    m = np.zero(n,n)
    unity_vert = 0.40 / n
    unity_hor  = 1 / n

    for i in range(0,n):
        for j in range(0, n):
            y = unity_vert * i
            x = unity_hor * j
            lam = find_line(x, y, p)
            length = length_line(lam)
            m[i,j] = 1.0 + 1/2 * length**2

    return m


def find_line(x, y, p) :
    lam = 0
    
    while f(x) < y :
        #mise a jour de f(x) fonction Romain
        lam = lam + p

    return lam

