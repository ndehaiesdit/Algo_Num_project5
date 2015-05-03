import numpy as np
import re
import matplotlib.pyplot as mp
import interpolation as ip

#-----------------------------------#
# Airfoil : load profile of a wing  #
#-----------------------------------#

def load_foil(file):
    f = open(file, 'r')
    matchline = lambda line: re.match(r"\s*([\d\.-]+)\s*([\d\.-]+)", line)
    extra  = [];    intra = []
    rextra = False; rintra = False
    for line in f:
        m = matchline(line)
        if (m != None) and not(rextra):
            rextra = True
        if (m != None) and rextra and not(rintra):
            extra.append(m.groups())
        if (m != None) and rextra and rintra:
            intra.append(m.groups())
        if (m == None) and rextra:
            rintra = True
    ex = np.array(map(lambda t: float(t[0]),extra))
    ey = np.array(map(lambda t: float(t[1]),extra))
    ix = np.array(map(lambda t: float(t[0]),intra))
    iy = np.array(map(lambda t: float(t[1]),intra))
    return(ex,ey,ix,iy)


def point_to_function(ax,ay,res):
    """Computes the polynomials coefficients corresponding with the value of its second derivate"""
    n = len(ax)
    Pol = np.zeros([n-1,4])
    for i in np.arange(0,n-1,1):
        tmp1 = (res[i+1]-res[i])/(ax[i+1]-ax[i])
        tmp2 = (ay[i+1]-ay[i])/(ax[i+1]-ax[i])
        Pol[i][0] = (1./6.)*tmp1 
        Pol[i][1] = (1./2.)*(res[i]-ax[i]*tmp1)
        Pol[i][2] = tmp2 - Pol[i][0]*(ax[i+1]**3-ax[i]**3)/(ax[i+1]-ax[i]) - Pol[i][1]*(ax[i+1]**2-ax[i]**2)/(ax[i+1]-ax[i])
        Pol[i][3] = ay[i] - Pol[i][0]*ax[i]**3 - Pol[i][1]*ax[i]**2 - Pol[i][2]*ax[i]

    return Pol

def general_function(x,ax,Pol):
    """Piecewise function of polynomials"""
    n = len(ax)
    for i in np.arange(0,n-1,1):
        if(x <= ax[i+1]):
            return Pol[i][0]*(x**3)+Pol[i][1]*(x**2)+Pol[i][2]*x+Pol[i][3]
    return Pol[n-1][0]*(x**3)+Pol[n-1][1]*(x**2)+Pol[n-1][2]*x+Pol[n-1][3]



(ex,ey,ix,iy) = load_foil("boe103.dat")
n = len(ex)
x = np.arange(0., 1.01, 0.01)

resulte = ip.second_derivate(ex,ey,0,1e30,n)
Pole = point_to_function(ex,ey,resulte)

resulti = ip.second_derivate(ix,iy,0,1e30,n)
Poli = point_to_function(ix,iy,resulti)
