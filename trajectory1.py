# coding:utf-8

import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import electro_module as pde

class CartesianGrid:                                                            #基本構造

    def __init__(self, nx=150, ny=150, nz=150, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=-5, zmax=25):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.ntotal = nx*ny*nz

        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

        self.dx = (xmax - xmin)/(nx - 1)
        self.dy = (ymax - ymin)/(ny - 1)
        self.dz = (zmax - zmin)/(nz - 1)

        self.x = np.arange(xmin, xmax + 0.5*self.dx, self.dx)
        self.y = np.arange(ymin, ymax + 0.5*self.dy, self.dy)
        self.z = np.arange(zmin, zmax + 0.5*self.dz, self.dz)

    def create_field(self):
        return np.zeros((self.nx, self.ny), dtype=np.float)                     #条件を入れるための格納庫、配列

    def create_meshgrid(self):                                                  #軸設定、max, minで表示する範囲を決定
        return np.meshgrid(self.x, self.y, self.z)

with open('electrode.binaryfile', 'rb') as lens:
    V = pickle.load(lens)

mesh = CartesianGrid()
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
I = I*(1/13)

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

def calc_trajectory(itera):
    

    for a in range(itera):
        data = [[],[],[]]
        zlist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
        ylist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
        vzlist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]

        Az = Ez*(q/m)
        Ay = Ey*(q/m)
        z = np.zeros(13)
        y = np.linspace(-4.5, 4.5, 13)

        for b in range(13):                         #軌道計算
            
            t = 0
            vz0 = 0
            vy0 = 0
            z0 = z[b]
            y0 = y[b]

            vzlist[b].append(vz0)
            zlist[b].append(z0)
            ylist[b].append(y0)

            while mesh.zmin+1<=z0<=mesh.zmax-0.001 and -11.5<=y0<=11.5:
                
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

            data[0].append(vz0)
            data[1].append(vy0)
            data[2].append(y0)

        for i in range(mesh.nx -1):                                 #電場リセット
            for j in range(mesh.nz -1):
                Ey[i, j] = -(V[i, j] - V[i+1, j])/delta_y
                Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

        iteration = 0
        for i in range(int((0-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin)), int((12-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin)-1)):  #軌道から電場
            alpha = 0.5 + iteration*0.5/int(12*mesh.nz/(mesh.zmax-mesh.zmin))
            iteration = iteration + 1
            
            y0list = []
            vz0list = []

            for j in range(13):
                y0list.append(getNearestValue(zlist, ylist, i*mesh.dz, j))
                vz0list.append(getNearestValue(zlist, vzlist, i*mesh.dz, j))
            
            r = np.max(y0list)
            v = np.mean(vz0list)

            for k in range(int((mesh.ymax-12)*mesh.ny/(mesh.ymax-mesh.ymin))+1, int((mesh.ymax+12)*mesh.ny/(mesh.ymax-mesh.ymin))):
                if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= k <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)): 
                    Ey[k-1, i] = Ey[k-1, i] + alpha*SpaceChargeEffect1(I, mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
                else:
                    Ey[k-1, i] = Ey[k-1, i] + alpha*SpaceChargeEffect2(I, mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))

        for i in range(int((12-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin)), int((mesh.zmax-mesh.zmin)*mesh.nz/(mesh.zmax-mesh.zmin))-1):  #軌道から電場
            y0list = []
            vz0list = []
            for j in range(11):
                y0list.append(getNearestValue(zlist, ylist, i*mesh.dz, j))
                vz0list.append(getNearestValue(zlist, vzlist, i*mesh.dz, j))
            
            r = np.max(y0list)
            v = np.mean(vz0list)

            for k in range(int((mesh.ymax-12)*mesh.ny/(mesh.ymax-mesh.ymin)), int((mesh.ymax+12)*mesh.ny/(mesh.ymax-mesh.ymin))):
                if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= k <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)): 
                    Ey[k, i] = Ey[k, i] + SpaceChargeEffect1(I, mesh.dy, r, v, k+1-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
                else:
                    Ey[k, i] = Ey[k, i] + SpaceChargeEffect2(I, mesh.dy, r, v, k+1-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))          

    with open('particles1.binaryfile', 'wb') as particles:
        pickle.dump(data, particles)
    
    return zlist, ylist
