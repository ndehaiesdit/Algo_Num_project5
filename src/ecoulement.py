import numpy as np
import matplotlib.pyplot as mp
import numpy.random as npr
import interpolation as ip
import integration as ig

# Tracement de l'aile

xhaut = np.array([0., 0.0177, 0.0532, 0.0952, 0.1419, 0.1903, 0.2387, 0.5838, 0.6435, 0.7112, 0.7709, 0.8306, 0.8854, 0.9532, 0.9870])

yhaut = np.array([0., 0.0274, 0.0440, 0.0614, 0.0730, 0.0797, 0.0847, 0.0779, 0.0730, 0.0647, 0.0548, 0.0456, 0.0348, 0.0191, 0.0116])

y2haut = ip.second_derivate(xhaut, yhaut, 0, 0, 999999)

n = np.size(xhaut)

y_interhaut = np.zeros(n)

for i in range(0, n - 1):
    y_interhaut[i] = ip.cubic_spline_meth(xhaut, yhaut, y2haut, xhaut[i])

    

xbas = np.array([0., 0.0096, 0.0419, 0.0951, 0.1370, 0.2016, 0.2629, 0.3129, 0.3838, 0.4370, 0.4870, 0.5887, 0.6693, 0.7351, 0.8064, 0.8709,
                  0.9241, 0.9693, 0.9903])

ybas = np.array([0., -0.0099, -0.0116, -0.0074, -0.0033, 0.0041, 0.0099, 0.0141, 0.0182, 0.0199, 0.0207, 0.0215, 0.0199, 0.0157, 0.0116, 0.0066,
                  0., -0.0041, -0.0074])

y2bas = ip.second_derivate(xbas, ybas, 0, 0, 999999)

n = np.size(xbas)

y_interbas = np.zeros(n)

for i in range(0, n - 1):
    y_interbas[i] = ip.cubic_spline_meth(xbas, ybas, y2bas, xbas[i])

mp.plot(xhaut,y_interhaut, 'g')
mp.plot(xbas, y_interbas, 'r')
mp.axis([0., 0.95, -0.10, 0.30])
mp.show()    


#Tracement des lignes de courants

def ligne(x, y, lamb, h):
    n = np.size(x)
    y_ligne = np.zeros(n)

    for i in range(0, n-1):
        y_ligne[i] = (1 - lamb) * y[i] + lamb*3*h

    return y_ligne


mp.plot(xhaut,y_interhaut, 'g')

lamb = 0.
while lamb < 1.:
    lamb = lamb + 0.03
    mp.plot(xhaut, ligne(xhaut, y_interhaut, lamb, 0.0872), 'black')
    
mp.plot(xbas, y_interbas, 'r')

lamb = 0.
while lamb < 1.:
    lamb = lamb + 0.08
    mp.plot(xbas, ligne(xbas, y_interbas, lamb, -0.0116), 'black')
    
    mp.axis([0., 0.95, -0.10, 0.30])
mp.show()          
