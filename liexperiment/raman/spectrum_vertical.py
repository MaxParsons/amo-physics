'''
Created on Jan 19, 2015

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
import liexperiment.traps
import amo.core.physicalconstants
import amo.quantum.trapstates
c = amo.core.physicalconstants.PhysicalConstantsSI
li = amo.core.physicalconstants.LithiumSixSI
from liexperiment.traps.lattice import composite_vertical_lattice

lat = composite_vertical_lattice()
lat.power_NE = 13
lat.power_NW = 21
lat.power_Ve = 10
zmin = -12.3e-6
zmax = -10.7e-6
npoints = 1000
zs = np.linspace(zmin, zmax, npoints)
depths = lat.axial_depth(zs)
n_states = 50

trap = amo.quantum.trapstates.Trap1D(zs, depths, li.mass, n_states)
trap.plot_eigenstates(skip=False, energy_rescale=1.0 / c.h / 1.0e6, x_rescale=1.0e6)
plt.show()
