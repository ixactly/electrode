# coding:utf-8

import numpy as np
import electro_module as pde
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm

with open('electrode.binaryfile', 'rb') as lens:
    V = pickle.load(lens)

mesh = pde.CartesianGrid()
V_pre = V[:, int(mesh.ny/2), :]
V = V_pre.transpose()

m = (40e-3)/(6.0*10e+23)
q = 1.6e-19
H = 1e-10

Ez = np.empty((mesh.nz-1, mesh.nx-1))
Ey = np.empty((mesh.nz-1, mesh.nx-1))

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

for i in range(mesh.nx -1):
    for j in range(mesh.nz -1):
        Ey[i, j] = (V[i+1, j] - V[i, j])/delta_y
        Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

Az = Ez*(q/m)
Ay = Ey*(q/m)

def Runge_Kutta(x0, a, v, h):

    k1 = v
    k2 = v+a*h/2
    k3 = v+a*h/2
    k4 = v+a*h

    x = x0 + 1000*(k1+2*k2+2*k3+k4)*h/6
    return x

def getNearestValue(zlist, ylist, value):                 #z方向の値を返して，対応するy方向の値に代入

    idx = np.abs(np.asarray(zlist) - value).argmin()
    number = zlist.index(idx)

    return ylist[number]

def SpaceChargeEffect(I, dy, y0, r, y0list, i, j):
        z = np.sqrt(r**2-(ylist[i])**2)
        d = y0 - y0list[i]
        pi = np.pi
        eps = 8.85418782e-12

        sigma = I/(4*pi*(r*10^-3)**2)

        E = dy*sigma*z/(2*pi*eps*d*np.sqrt(d**2+z**2))
        return E

zlist = [[],[],[],[],[],[],[],[],[],[],[]]
ylist = [[],[],[],[],[],[],[],[],[],[],[]]

z = np.zeros(11)
y = np.linspace(-4.6, 4.6, 11)

for i in range(11):
    t = 0
    vz0 = 0
    vy0 = 0
    z0 = z[i]
    y0 = y[i]

    while mesh.zmin+1<=z0<=mesh.zmax-1 and mesh.ymin+1<=y0<=mesh.ymax-1:
        t += H

        az = Az[int((mesh.ny-1)*(mesh.ymax-y0)/(mesh.ymax-mesh.ymin)), int((mesh.nz-1)*(z0-mesh.zmin)/(mesh.zmax - mesh.zmin))]
        ay = Ay[int((mesh.ny-1)*(mesh.ymax-y0)/(mesh.ymax-mesh.ymin)), int((mesh.nz-1)*(z0-mesh.zmin)/(mesh.zmax - mesh.zmin))]

        vz0 += az*H
        vy0 += ay*H

        z0 = Runge_Kutta(z0, az, vz0, H)
        y0 = Runge_Kutta(y0, ay, vy0, H)

        zlist[i].append(z0)
        ylist[i].append(y0)
        
zform1 = [mesh.zmax, 0, 0, 2, 2, mesh.zmax] 
yform1 = [20, 20, 5, 5, 18.5, 18.5]

zform2 = [mesh.zmax, 10, 10, 12, 12, mesh.zmax]
yform2 = [13, 13, 5, 5, 11.5, 11.5]

yreform1 = [-20, -20, -5, -5, -18.5, -18.5]
yreform2 = [-13, -13, -5, -5, -11.5, -11.5]

fig = plt.figure()

for i in range(11):
        plt.plot(zlist[i], ylist[i], color="r")

plt.plot(zform1, yform1, color="k")
plt.plot(zform2, yform2, color="k")
plt.plot(zform1, yreform1, color="k")
plt.plot(zform2, yreform2, color="k")

plt.xlim([mesh.zmin, mesh.zmax])
plt.ylim([mesh.xmin, mesh.xmax])

plt.show()
print(z, y)