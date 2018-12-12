# coding:utf-8

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import csc_matrix
import time


class CartesianGrid:                                                            #基本構造
    """
        Simple class to generate a computational grid and apply boundary conditions
    """

    def __init__(self, nx=150, ny=150, nz=150, xmin=-25, xmax=25, ymin=-25, ymax=25, zmin=-50, zmax=50):
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

    def set_boundary_condition1(self, V_1, r, eps=0.1):                                  #円柱状にデータを配列する
        a = r*r - V_1                                                             #rをmeshgridに合わせて考えること
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2 - a
        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] > V_1 + eps or Z[i, j] < V_1 - eps:
                    Z[i, j] = 0

        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def set_boundary_condition2(self, phi, V_0=1e-100):
        for i in range(self.nx):
            for j in range(self.ny):
                if phi[0, i, j] != 0.0:            #修正
                    phi[0, i, j] = V_0
        return phi

    def make_einzel_lens(self, Z_V, Z_0, z1, z2, z3, z4, z5, z6):
        Zeros = np.zeros((1, self.nx, self.ny))
        Z_new = np.zeros((1, self.nx, self.ny))
        for i in range(int(self.nz*z1/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Zeros))

        for i in range(int(self.nz*z1/(self.zmax - self.zmin)), int(self.nz*z2/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Z_0))

        for i in range(int(self.nz*z2/(self.zmax - self.zmin)), int(self.nz*z3/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Zeros))

        for i in range(int(self.nz*z3/(self.zmax - self.zmin)), int(self.nz*z4/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Z_V))

        for i in range(int(self.nz*z4/(self.zmax - self.zmin)), int(self.nz*z5/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Zeros))

        for i in range(int(self.nz*z5/(self.zmax - self.zmin)), int(self.nz*z6/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, Z_0))

        for i in range(int(self.nz*z6/(self.zmax - self.zmin)), int(self.nz)-1):
            Z_new = np.vstack((Z_new, Zeros))

        return Z_new

    def convert_to_1d_array(self, x):
        return x.reshape(self.ntotal, 1)

    def convert_to_3d_array(self, x):
        return x.reshape(self.nx, self.ny, self.nz)
