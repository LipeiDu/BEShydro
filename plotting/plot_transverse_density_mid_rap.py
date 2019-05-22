import numpy as np
import matplotlib.pyplot as plt
import sys

plt.style.use('classic')
filename = sys.argv[1]

nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])

dx = 0.1
dy = 0.1
dz = 0.5

xmax = dx*((nx-1)/2)
ymax = dy*((ny-1)/2)
x, y, z, v = np.loadtxt(filename, unpack=True)

xnew = np.linspace(-xmax, xmax, num=nx)
ynew = np.linspace(-ymax, ymax, num=ny)

vnew = v.reshape(nz, ny, nx).transpose()

front = vnew[:, :, 0].transpose()
back = vnew[:, :, -1].transpose()
middle = vnew[:, :, ((nz-1)/2)].transpose()

fig, ax = plt.subplots()
cset1 = ax.pcolormesh(xnew, ynew, middle)
plt.xlabel('x [fm]')
plt.ylabel('y [fm]')
plt.title('transverse energy density')
cbar = plt.colorbar(cset1)
cbar.ax.set_ylabel('energy density [fm^-4]')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
