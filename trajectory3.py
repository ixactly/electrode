# coding:utf-8

import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import three_d_pde as pde

with open('einzel_lens.binaryfile', 'rb') as lens:
    V = pickle.load(lens)

mesh = pde.CartesianGrid()
V_pre = V[:, int(mesh.ny/2), :]
V = V_pre.transpose()

m = (40e-3)/(6.0*1e+23)
q = 1.6e-19
H = 1e-10
eps = 8.85418782e-12
V_extract = 2000
d = 1
J = (4*eps*(V_extract**(3/2))*np.sqrt(2*q/m))/(9*d**2)
I = J*np.pi*d**2/4
I = I*0.65
itera = 4

Ez = np.empty((mesh.nz-1, mesh.nx-1))
Ey = np.empty((mesh.nz-1, mesh.nx-1))

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

with open('particles2.binaryfile', 'rb') as particle:
    data = pickle.load(particle)

m = (40e-3)/(6.0*1e+23)
q = 1.6e-19
H = 1e-10
eps = 8.85418782e-12
V_extract = 2000
d = 1
J = (4*eps*(V_extract**(3/2))*np.sqrt(2*q/m))/(9*d**2)
I = J*np.pi*d**2/4
I = I*0.65
itera = 20

Ez = np.empty((mesh.nz-1, mesh.nx-1))
Ey = np.empty((mesh.nz-1, mesh.nx-1))

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

for i in range(mesh.nx -1):                                 #電場リセット
    for j in range(mesh.nz -1):
        Ey[i, j] = -(V[i, j] - V[i+1, j])/delta_y
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

    sigma = I/(v*pi*(r*1e-3)**2)

    E = -sigma*i*dy*1e-3/(2*eps)
    return E

def SpaceChargeEffect2(I, dy, r, v, i):
    pi = np.pi
    eps = 8.85418782e-12

    sigma = I/v

    E = -sigma/(2*pi*eps*(i*dy*1e-3))
    return E
def calc_trajectory(itera, data):
    data1 = [[],[],[]]

    for a in range(itera):
        zlist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
        ylist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
        vzlist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]

        Az = Ez*(q/m)
        Ay = Ey*(q/m)

        for b in range(13):                         #軌道計算
            
            t = 0
            vz0 = data[0][b]
            vy0 = data[1][b]
            z0 = mesh.zmin
            y0 = data[2][b]

            while mesh.zmin-0.3<=z0<=mesh.zmax and -18.5<=y0<=18.5:
                
                t += H

                az = Az[int((mesh.ny-1)*(mesh.ymax-y0)/(mesh.ymax-mesh.ymin)), int((mesh.nz-1)*(z0-mesh.zmin)/(mesh.zmax - mesh.zmin))]
                ay = Ay[int((mesh.ny-1)*(mesh.ymax-y0)/(mesh.ymax-mesh.ymin)), int((mesh.nz-1)*(z0-mesh.zmin)/(mesh.zmax - mesh.zmin))]

                vz0 += az*H
                vy0 += ay*H

                z0 = Runge_Kutta(z0, az, vz0, H)
                y0 = Runge_Kutta(y0, ay, vy0, H)

                vzlist[b].append(vz0)
                zlist[b].append(z0)
                ylist[b].append(y0)

            data1[0].append(vz0)
            data1[1].append(vy0)
            data1[2].append(y0)

        for i in range(mesh.nx -1):                                 #電場リセット
            for j in range(mesh.nz -1):
                Ey[i, j] = -(V[i, j] - V[i+1, j])/delta_y
                Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

        for i in range(int(mesh.nz-1)):  #軌道から電場
            co = 4
            y0list = []
            vz0list = []
            for j in range(13):
                if zlist[j] != []:
                    y0list.append(getNearestValue(zlist, ylist, mesh.zmin+i*mesh.dz, j))
                    vz0list.append(getNearestValue(zlist, vzlist, mesh.zmin+i*mesh.dz, j))
            
            r = np.max(y0list)
            v = np.mean(vz0list)

            for j in range(2,11):
                if y0list[j] > 18.5 or y0list[j] < -18.5:
                    co += 1
                
            for k in range(int((mesh.ymax-12)*mesh.ny/(mesh.ymax-mesh.ymin))+1, int((mesh.ymax+12)*mesh.ny/(mesh.ymax-mesh.ymin))):
                if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= k <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)): 
                    Ey[k-1, i] = Ey[k-1, i] + SpaceChargeEffect1(I*(1-co/13), mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
                else:
                    Ey[k-1, i] = Ey[k-1, i] + SpaceChargeEffect2(I*(1-co/13), mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
            
    zform1 = [-45, -25, -25, -45, -45]
    zform2 = [-20, 10, 10, -20, -20]
    zform3 = [15, 45, 45, 15, 15]

    yform1 = [-20, -20, -18.5, -18.5, -20]
    yform2 = [20, 20, 18.5, 18.5, 20]

    for i in range(2, 11):
        plt.plot(zlist[i], ylist[i], color="r")

    plt.plot(zform1, yform1, color="k")
    plt.plot(zform1, yform2, color="k")
    plt.plot(zform2, yform1, color="k")
    plt.plot(zform2, yform2, color="k")
    plt.plot(zform3, yform1, color="k")
    plt.plot(zform3, yform2, color="k")

    plt.xlim([mesh.zmin, mesh.zmax])
    plt.ylim([mesh.xmin, mesh.xmax])

    plt.show()

    with open('particles3.binaryfile', 'wb') as particles:
        pickle.dump(data, particles)
calc_trajectory(itera, data)