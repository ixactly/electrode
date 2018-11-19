# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
import practice as pde
import pickle
from matplotlib import cm

with open('electrode.binaryfile', 'rb') as electrode:
    V = pickle.load(electrode)

mesh = pde.CartesianGrid()
V_pre = V[:, int(mesh.ny/2), :]
two_d_V = V_pre.transpose()

x, y = np.meshgrid(mesh.z, mesh.x)
fig, ax = plt.subplots()
surf = ax.contourf(x, y, two_d_V, cmap=cm.coolwarm)

fig.colorbar(surf, shrink=0.5, aspect=5).set_label("potential[V]")
plt.gca().set_aspect('equal', adjustable='box')

plt.show()
