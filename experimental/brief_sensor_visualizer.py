import numpy as np
from scipy.stats import mstats
import scipy.misc
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("Agg")
from matplotlib.ticker import MultipleLocator
import matplotlib.animation as manimation
import sys, traceback, os, math, argparse

class CONF: pass
conf = CONF()
sensors = np.loadtxt("../output/sensors_Nx1552_Ny1328.txt")
conf.subdomain_y = 16
conf.subdomain_x = 16
Nx = 97
Ny = 83

print("number of sensors = " + str(sensors.shape[0]))
print("number of subdomains = " + str(Nx * Ny))
print("number of nodal point = " +
        str(Nx * Ny * conf.subdomain_x * conf.subdomain_y))

fig = plt.figure()
ax = fig.gca()
ax.scatter(sensors[:,0], sensors[:,1])
ax.grid(which="minor", color="r", linestyle=":")
ax.xaxis.set_minor_locator(MultipleLocator(conf.subdomain_x))
ax.yaxis.set_minor_locator(MultipleLocator(conf.subdomain_y))
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim((0, Nx * conf.subdomain_x))
ax.set_ylim((0, Ny * conf.subdomain_y))
plt.title("Sensor locations")
plt.tight_layout()
plt.show()
plt.savefig("../output/1.jpg")


