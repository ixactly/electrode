# coding:utf-8

import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import electro_module as pde

with open('electrode.binaryfile', 'rb') as lens:
    V = pickle.load(lens)

mesh = pde.CartesianGrid()
V_pre = V[:, int(mesh.ny/2), :]
V = V_pre.transpose()

m = (40e-3)/(6.0*10e+23)
q = 1.6e-19
H = 1e-10
I = 1

Ez = np.empty((mesh.nz-1, mesh.nx-1))
Ey = np.empty((mesh.nz-1, mesh.nx-1))

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

for i in range(mesh.nx -1):
    for j in range(mesh.nz -1):
        Ey[i, j] = (V[i+1, j] - V[i, j])/delta_y
        Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

def Runge_Kutta(x0, a, v, h):

    k1 = v
    k2 = v+a*h/2
    k3 = v+a*h/2
    k4 = v+a*h

    x = x0 + 1000*(k1+2*k2+2*k3+k4)*h/6
    return x

def getNearestValue(zlist, ylist, value, i):                 #z方向の値を返して，対応するy方向の値に代入

    idx = np.abs(np.asarray(zlist[i]) - value).argmin()
    return ylist[i][idx]

def SpaceChargeEffect1(I, dy, r, v, i):
    pi = np.pi
    eps = 8.85418782e-12

    sigma = I/(v*pi*(r*(1e-3))**2)

    E = sigma*i*dy*1e-3/2*eps
    return E

def SpaceChargeEffect2(I, r, v):
    pi = np.pi
    eps = 8.85418782e-12

    sigma = I/(v*pi*(r*1e-3)**2)

    E = sigma*r/2*eps
    return E
    
zlist = [[],[],[],[],[],[],[],[],[],[],[]]
ylist = [[],[],[],[],[],[],[],[],[],[],[]]
vzlist = [[],[],[],[],[],[],[],[],[],[],[]]

z = np.zeros(11)
y = np.linspace(-4.6, 4.6, 11)

# main code
for a in range(2):
    Az = Ez*(q/m)
    Ay = Ey*(q/m)
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

            vzlist[i].append(vz0)
            zlist[i].append(z0)
            ylist[i].append(y0)

    for i in range(int((0-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin)), int((mesh.zmax-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin))): 
        for j in range(11):
            y0list = []
            y0list.append(getNearestValue(zlist[j], ylist, i*mesh.dz, j))
            vz0list = []
            vz0list.append(getNearestValue(zlist[j], vzlist, i*mesh.dz, j))

        r = np.max(y0list)
        E = 0
        v = np.mean(vz0list)

        for k in range(int((mesh.ymax-12)*mesh.ny/(mesh.ymax-mesh.ymin)), int((mesh.ymax+12)*mesh.ny/(mesh.ymax-mesh.ymin))+1):
            if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= k <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)): 
                Ey[i, k] = Ey[i, k] + SpaceChargeEffect1(I, mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))

            else:
                Ey[i, k] = Ey[i, k] + SpaceChargeEffect2(I, r, v)

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
