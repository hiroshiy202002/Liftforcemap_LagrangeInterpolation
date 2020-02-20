#! /usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME = sys.argv[1]
N = 160
NODES = N+1

data = np.loadtxt(FILENAME, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
F = data[2].reshape(NODES, NODES)

fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(1, 1, 1)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
### contour
cont = ax.contour(Y, Z, F, levels=[0.0], colors="#000000", linewidths=0.4)
cont.clabel(fmt="%1.1f", fontsize=9)
### color map
val1 = np.amin(F)
val2 = np.amax(F)
c = ax.pcolormesh(Y, Z, F, cmap="jet", vmin=val1, vmax=val2)
plt.colorbar(c, format="%1.0e")
ax.set_aspect("equal")

plt.tight_layout()
plt.suptitle(FILENAME, fontsize=20)
plt.subplots_adjust(top=0.94)

plt.savefig("output.ps")
#   plt.show()
