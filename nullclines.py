#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME1 = "IntpltdFr.dat"
FILENAME2 = "IntpltdFt.dat"
FILENAME3 = "EquilibriumPoints.dat"
N = 160
NODES = N+1

data = np.loadtxt(FILENAME1, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
Fr = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME2, unpack=True)
Ft = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME3, unpack=True)
Ye = data[0]
Ze = data[1]

fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(1, 1, 1)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
### contour
cont = ax.contour(Y, Z, Fr, levels=[0.0], colors="#000000", linestyles="dashed")
cont = ax.contour(Y, Z, Ft, levels=[0.0], colors="#000000")
### points
plt.scatter(Ye, Ze, marker="x", c="red")
ax.set_aspect("equal")

plt.tight_layout()
plt.suptitle("Nullclines", fontsize=20)
plt.subplots_adjust(top=0.94)

plt.savefig("output.ps")
#   plt.show()
