# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
import time
import pickle

class CartesianGrid:                                                            #基本構造

    def __init__(self, nx=300, ny=300, nz=300, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=25, zmax=75):
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


with open('particles.binaryfile', 'rb') as particle:
    data = pickle.load(particle)

mesh = CartesianGrid()

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

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

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

    Ez = np.zeros((mesh.nz-1, mesh.nx-1))
    Ey = np.zeros((mesh.nz-1, mesh.nx-1))

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
            z0 = 25
            y0 = data[2][b]

            while mesh.zmin-0.3<=z0<=mesh.zmax-0.1 and -11.5<=y0<=11.5:
                
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

        Ez = np.zeros((mesh.nz-1, mesh.nx-1))
        Ey = np.zeros((mesh.nz-1, mesh.nx-1))

        for i in range(int(mesh.nz-1)):  #軌道から電場
            co = 0
            y0list = []
            vz0list = []
            for j in range(13):
                y0list.append(getNearestValue(zlist, ylist, mesh.zmin+i*mesh.dz, j))
                vz0list.append(getNearestValue(zlist, vzlist, mesh.zmin+i*mesh.dz, j))
            
            r = np.max(y0list)
            v = np.mean(vz0list)

            for j in range(13):
                if y0list[j] > 11.2 or y0list[j] < -11.2:
                    co += 1
                print(y0list[j])
            print(co)
            print(mesh.zmin+mesh.dz*i)
            for k in range(int((mesh.ymax-12)*mesh.ny/(mesh.ymax-mesh.ymin))+1, int((mesh.ymax+12)*mesh.ny/(mesh.ymax-mesh.ymin))):
                if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= k <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)): 
                    Ey[k-1, i] = Ey[k-1, i] + SpaceChargeEffect1(I*(1-co/13), mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
                else:
                    Ey[k-1, i] = Ey[k-1, i] + SpaceChargeEffect2(I*(1-co/13), mesh.dy, r, v, k-int(mesh.ymax*mesh.ny/(mesh.ymax-mesh.ymin)))
            
    zform1 = [mesh.zmax+40, 0, 0, 2, 2, mesh.zmax]
    yform1 = [20, 20, 5, 5, 18.5, 18.5]

    zform2 = [mesh.zmax+40, 10, 10, 12, 12, mesh.zmax]
    yform2 = [13, 13, 5, 5, 11.5, 11.5]

    yreform1 = [-20, -20, -5, -5, -18.5, -18.5]
    yreform2 = [-13, -13, -5, -5, -11.5, -11.5]

    for i in range(13):
        plt.plot(zlist[i], ylist[i], color="r")

    plt.plot(zform1, yform1, color="k")
    plt.plot(zform2, yform2, color="k")
    plt.plot(zform1, yreform1, color="k")
    plt.plot(zform2, yreform2, color="k")

    plt.xlim([mesh.zmin, mesh.zmax])
    plt.ylim([mesh.xmin, mesh.xmax])

    plt.show()

    with open('particles2.binaryfile', 'wb') as particles:
        pickle.dump(data1, particles)
        
calc_trajectory(itera, data)