# coding utf-8

import numpy as np
import math
from random import random
import matplotlib.pyplot as plt

m = (40e-3)/(6.0*1e+23)
q = 1.6e-19
H = 1e-10
eps = 8.85418782e-12
V0 = 2000
d = 1
J = (4*eps*(V0**(3/2))*np.sqrt(2*q/m))/(9*d**2)
I = J*np.pi*d**2/4
I = I/2
K = I*np.sqrt(m/q)/(np.pi*eps*(2*V0)**(3/2))


def RungeKutta(z0, u0, r0, h, K):
    #dr/dz = u, du/dz = K/2r → dr/dz = f(z, u, r)  du/dz = g(z, u, r)
    liist = []
    au = K/(2*r0)
    ar = u0

    bu = K/(2*(r0+h/2))
    br = u0 + au*h/2

    cu = K/(2*(r0+h/2))
    cr = u0 + bu*h/2
    
    du = K/(2*(r0+h))
    dr = u0 + cu*h

    r1 = r0 + (ar + br + cr + dr)*h/6
    u1 = u0 + (au + bu + cu + du)*h/6

    liist.append(u1)
    liist.append(r1)

    return liist

z0 = 0
u0 = 0
r0 = d/2
h = 0.01
zlist1 = []
rlist1 = []
zlist1.append(z0)
rlist1.append(r0)

while z0 < 50:
    z0 = z0 + h
    data = RungeKutta(z0, u0, r0, h, K)
    u0 = data[0]
    r0 = data[1]

    zlist1.append(z0)
    rlist1.append(r0)
    
V0 = 10000
d = 2
J = (4*eps*(V0**(3/2))*np.sqrt(2*q/m))/(9*d**2)
I = J*np.pi*d**2/4 
I = 0.01
K = I*np.sqrt(m/q)/(np.pi*eps*(2*V0)**(3/2))
z0 = 0
u0 = 0
r0 = d/2
h = 0.01
zlist2 = []
rlist2 = []
zlist2.append(z0)
rlist2.append(r0)

while z0 < 50:
    z0 = z0 + h
    data = RungeKutta(z0, u0, r0, h, K)
    u0 = data[0]
    r0 = data[1]

    zlist2.append(z0)
    rlist2.append(r0)
plt.title("V=2[kV], I=0.1[mA]       V=10[kV], I=10[mA]")
plt.plot(zlist1, rlist1)
plt.plot(zlist2, rlist2)
plt.show()
    