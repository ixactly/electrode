# coding:utf-8

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import csc_matrix
import time
import pickle

class CartesianGrid:                                                            #基本構造
    """
        Simple class to generate a computational grid and apply boundary conditions
    """

    def __init__(self, nx=150, ny=150, nz=150, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=-5, zmax=35):
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

    def set_boundary_circle(self, V_1, r1, r2):                                  #円柱状にデータを配列する                                                             #rをmeshgridに合わせて考えること
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r1**2:
                    Z[i, j] = 0
                if Z[i, j] >= r2**2:
                    Z[i, j] = 0
                if Z[i, j] != 0:
                    Z[i, j] = V_1

        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def set_boundary_2circle(self, V_1, V_2, r1, r2, r3, r4):
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r1**2:
                    Z[i, j] = 0
                if r2**2 <= Z[i, j] <= r3**2:
                    Z[i, j] = 0
                if Z[i, j] >= r4**2:
                    Z[i, j] = 0
                if r1**2 <= Z[i, j] <= r2**2:
                    Z[i, j] = V_1
                if r3**2 <= Z[i, j] <= r4**2:
                    Z[i, j] = V_2
        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def make_cylinder(self, a, b, c, d, e, z1, z2, z3, z4):
        Z_new = e
        for i in range(int(self.nz*(z1-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, b))

        for i in range(int(self.nz*(z1-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z2-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, a))

        for i in range(int(self.nz*(z2-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z3-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, b))

        for i in range(int(self.nz*(z3-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z4-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, c))

        for i in range(int(self.nz*(z4-self.zmin)/(self.zmax - self.zmin)), int(self.nz)-1):
            Z_new = np.vstack((Z_new, d))

        return Z_new

    def convert_to_1d_array(self, x):
        return x.reshape(self.ntotal, 1)

    def convert_to_3d_array(self, x):
        return x.reshape(self.nx, self.ny, self.nz)

class coordinate:
    def __init__(self, i, j, k, mesh):
        self.p = i*mesh.nx*mesh.ny + j*mesh.nx + k
        self.ip1 = (i+1)*mesh.nx*mesh.ny + j*mesh.nx + k
        self.im1 = (i-1)*mesh.nx*mesh.ny + j*mesh.nx + k
        self.jp1 = i*mesh.nx*mesh.ny + (j+1)*mesh.nx + k
        self.jm1 = i*mesh.nx*mesh.ny + (j-1)*mesh.nx + k
        self.kp1 = i*mesh.nx*mesh.ny + j*mesh.nx + (k+1)
        self.km1 = i*mesh.nx*mesh.ny + j*mesh.nx + (k-1)

def calc_jacobi_matrix(mesh, Z):
    """
        Create sparse matrix for Jacobi method
    """

    A = lil_matrix((mesh.ntotal, mesh.ntotal))
    i, j, k = 0, 0, 0
    for i in range(1, mesh.nz-1):
        for j in range(1, mesh.ny-1):
            for k in range(1, mesh.nx-1):
                p = coordinate(i, j, k, mesh)

                if Z[i, j, k] != 0:                                             #0の時はe-10などで近似して、なんとかする
                    A[p.p, p.p] = 1.0

                else:
                    A[p.p, p.ip1] = 1/6
                    A[p.p, p.im1] = 1/6
                    A[p.p, p.jp1] = 1/6
                    A[p.p, p.jm1] = 1/6
                    A[p.p, p.kp1] = 1/6
                    A[p.p, p.km1] = 1/6

    i, j, k = 0, 0, 0
    p = coordinate(i, j, k, mesh)
    A[p.p, p.ip1] = 1/6
    A[p.p, p.jp1] = 1/6
    A[p.p, p.kp1] = 1/6

    i, j, k = mesh.nz-1, 0, 0
    p = coordinate(i, j, k, mesh)
    A[p.p, p.im1] = 1/6
    A[p.p, p.jp1] = 1/6
    A[p.p, p.kp1] = 1/6

    i, j, k = 0, mesh.ny-1, 0
    p = coordinate(i, j, k, mesh)
    A[p.p, p.ip1] = 1/6
    A[p.p, p.jm1] = 1/6
    A[p.p, p.kp1] = 1/6

    i, j, k = mesh.nz-1, mesh.ny-1, 0
    p = coordinate(i, j, k, mesh)
    A[p.p, p.im1] = 1/6
    A[p.p, p.jm1] = 1/6
    A[p.p, p.kp1] = 1/6

    i, j, k = 0, 0, mesh.nx-1
    p = coordinate(i, j, k, mesh)
    A[p.p, p.ip1] = 1/6
    A[p.p, p.jp1] = 1/6
    A[p.p, p.km1] = 1/6

    i, j, k = mesh.nz-1, 0, mesh.nx-1
    p = coordinate(i, j, k, mesh)
    A[p.p, p.im1] = 1/6
    A[p.p, p.jp1] = 1/6
    A[p.p, p.km1] = 1/6

    i, j, k = 0, mesh.ny-1, mesh.nx-1
    p = coordinate(i, j, k, mesh)
    A[p.p, p.ip1] = 1/6
    A[p.p, p.jm1] = 1/6
    A[p.p, p.km1] = 1/6

    i, j, k = mesh.nz-1, mesh.ny-1, mesh.nx-1
    p = coordinate(i, j, k, mesh)
    A[p.p, p.im1] = 1/6
    A[p.p, p.jm1] = 1/6
    A[p.p, p.km1] = 1/6

    # (Z,Y,X) = (i,0,0)
    i, j, k = 1, 0, 0
    for i in range(1, mesh.nz-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.kp1] = 1/6
    # (Z,Y,X) = (i,ymax,0)
    i, j, k = 1, mesh.ny-1, 0
    for i in range(1, mesh.nz-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.kp1] = 1/6
    # (Z,Y,X) = (i,0,zmax)
    i, j, k = 1, 0, mesh.nz-1
    for i in range(1, mesh.nz-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (i,ymax,zmax)
    i, j ,k = 1, mesh.ny-1, mesh.nz-1
    for i in range(1, mesh.nz-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.ip1] = 1/6
        A[p.p, p.im1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (0,j,0)
    i, j, k = 0, 1, 0
    for j in range(1, mesh.ny-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.kp1] = 1/6
    # (Z,Y,X) = (zmax,j,0)
    i, j, k = mesh.nz-1, 1, 0
    for j in range(1, mesh.ny-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.kp1] = 1/6
    # (Z,Y,X) = (0,j,xmax)
    i, j, k = 0, 1, mesh.nx-1
    for j in range(1, mesh.ny-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (zmax,j,zmax)
    i, j, k = mesh.nz-1, 1, mesh.ny-1
    for j in range(1, mesh.ny-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (0,0,k)
    i, j, k = 0, 0, 1
    for k in range(1, mesh.nx-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.kp1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (zmax,0,k)
    i, j, k = mesh.nz-1, 0, 1
    for k in range(1, mesh.nx-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.jp1] = 1/6
        A[p.p, p.kp1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (0,ymax,k)
    i, j, k = 0, mesh.ny-1, 1
    for k in range(1, mesh.nx-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.ip1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.kp1] = 1/6
        A[p.p, p.km1] = 1/6
    # (Z,Y,X) = (zmax,ymax,k)
    i, j, k = mesh.nz-1, mesh.ny-1, 1
    for k in range(1, mesh.nx-1):
        p = coordinate(i, j, k, mesh)
        A[p.p, p.im1] = 1/6
        A[p.p, p.jm1] = 1/6
        A[p.p, p.kp1] = 1/6
        A[p.p, p.km1] = 1/6

    for j in range(1, mesh.ny-1):
        for k in range(1, mesh.nx-1):
            i = 0
            p = coordinate(i, j, k, mesh)
            if Z[i, j, k] != 0:
                A[p.p, p.p] = 1.0
            else:
                A[p.p, p.ip1] = 1/6
                A[p.p, p.jp1] = 1/6
                A[p.p, p.jm1] = 1/6
                A[p.p, p.kp1] = 1/6
                A[p.p, p.km1] = 1/6

            i = mesh.nz-1
            p = coordinate(i, j, k, mesh)
            if Z[i, j, k] != 0:
                A[p.p, p.p] = 1.0
            else:
                A[p.p, p.im1] = 1/6
                A[p.p, p.jp1] = 1/6
                A[p.p, p.jm1] = 1/6
                A[p.p, p.kp1] = 1/6
                A[p.p, p.km1] = 1/6

    for k in range(1, mesh.nx-1):
        for i in range(1, mesh.nz-1):
            j = 1
            p = coordinate(i, j, k, mesh)
            A[p.p, p.ip1] = 1/6
            A[p.p, p.im1] = 1/6
            A[p.p, p.jp1] = 1/6
            A[p.p, p.kp1] = 1/6
            A[p.p, p.km1] = 1/6

            j = mesh.ny-1
            p = coordinate(i, j, k, mesh)
            A[p.p, p.ip1] = 1/6
            A[p.p, p.im1] = 1/6
            A[p.p, p.jm1] = 1/6
            A[p.p, p.kp1] = 1/6
            A[p.p, p.km1] = 1/6

    for i in range(1, mesh.nz-1):
        for j in range(1, mesh.ny-1):
            k = 0
            p = coordinate(i, j, k, mesh)
            A[p.p, p.ip1] = 1/6
            A[p.p, p.im1] = 1/6
            A[p.p, p.jp1] = 1/6
            A[p.p, p.jm1] = 1/6
            A[p.p, p.kp1] = 1/6

            k = mesh.nx-1
            p = coordinate(i, j, k, mesh)
            A[p.p, p.ip1] = 1/6
            A[p.p, p.im1] = 1/6
            A[p.p, p.jm1] = 1/6
            A[p.p, p.jp1] = 1/6
            A[p.p, p.km1] = 1/6

    return A.tocsc()

class IterationControl:
    """
        Class to control iteration loop
    """

    def __init__(self, max_iter, info_interval, tolerance):
        self.max_iter = max_iter
        self.info_interval = info_interval
        self.tolerance = tolerance
        self.eps = 1.0
        self.iter = 0

    def loop(self):
        self.iter += 1
        self.output_info()

        if self.eps < self.tolerance:
            return False
        elif self.iter > self.max_iter:
            print("max iteration reached")
            return False
        else:
            return True

    def calc_epsilon(self, dx):
        self.eps = np.max(abs(dx))

    def output_info(self):
        if self.iter % self.info_interval == 0:
            print("iter = %d, eps = %.3e" % (self.iter, self.eps))

#main code
def solve_eq():

    mesh = CartesianGrid()

    # extrode boundary condition
    A = mesh.set_boundary_circle(2000, 5, 19)
    B = mesh.set_boundary_circle(2000, 17.8, 19)
    C = mesh.set_boundary_2circle(1e-100, 2000, 5, 12.7, 17.8, 19)
    D = mesh.set_boundary_2circle(1e-100, 2000, 11.5, 12.7, 17.8, 19)
    E = mesh.set_boundary_circle(2000, 0, 19)
    electrode = mesh.make_cylinder(A, B, C, D, E, 0, 2, 12, 14)

    A = calc_jacobi_matrix(mesh, electrode)

    k = mesh.convert_to_1d_array(electrode)

    iter_control = IterationControl(25000, 100, 1e-3)

    while iter_control.loop():
        k_new = A.dot(k)
        iter_control.calc_epsilon(k_new - k)
        k, k_new = k_new, k

    # reshape for surface plotting
    k = mesh.convert_to_3d_array(k)
    return k

Potential = solve_eq()

with open('electrode.binaryfile', 'wb') as lens:
    pickle.dump(Potential, lens)
