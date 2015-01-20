'''
Created on Jan 19, 2015

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
import os.path
import amo.core.simulation
import liexperiment.traps
import amo.core.physicalconstants
import amo.quantum.trapstates
import amo.core.simulation
c = amo.core.physicalconstants.PhysicalConstantsSI
li = amo.core.physicalconstants.LithiumSixSI
from liexperiment.traps.lattice import composite_vertical_lattice

# lat = composite_vertical_lattice()
# lat.power_NE = 13
# lat.power_NW = 21
# lat.power_Ve = 10
# zmin = -12.3e-6
# zmax = -10.7e-6
# npoints = 1000
# zs = np.linspace(zmin, zmax, npoints)
# depths = lat.axial_depth(zs)
# n_states = 50
#
# trap = amo.quantum.trapstates.Trap1D(zs, depths, li.mass, n_states)
# trap.plot_eigenstates(skip=False, energy_rescale=1.0 / c.h / 1.0e6, x_rescale=1.0e6)
# plt.show()

class vertical_lattice_spectrum(amo.core.simulation.simulation):
    def shot(self, shotnumber, parameters):
        shot_info_string = str()
        for var in self.scan_variables:
            shot_info_string += "{0} = {1}\n".format(var, parameters[var])
        print shot_info_string
        lat = composite_vertical_lattice()
        lat.power_NE = parameters["power_NE"]
        lat.power_NW = parameters["power_NE"]
        lat.power_Ve = parameters["power_Ve"]
        zmin = parameters["zmin"]
        zmax = parameters["zmax"]
        npoints = parameters["npoints"]
        n_states = parameters["nstates"]
        zs = np.linspace(zmin, zmax, npoints)
        depths = lat.axial_depth(zs)
        trap = amo.quantum.trapstates.Trap1D(zs, depths, li.mass, n_states, npol=10)
        trap.compute_occupations();
        energies, spectrum = trap.simulate_spectrum(2.0 * np.pi / (0.671e-6))
        
        if self.config["is_debug"]:
            figstates, axstates = plt.subplots(1, 1)
            trap.plot_eigenstates(axstates, energy_rescale=1.0 / (c.h * 1.0e6), x_rescale=1.0e6)
            axstates.set_title("Trap Eigenstates\n" + shot_info_string)
            
        
        
dir = "C:\\Users\\Max\\Documents\\Dropbox\\Physics\\Lab\\Python Calculations\\VerticalLatticeSpectrum"        
parameters = {
              "power_NE" : np.linspace(0.0, 1, 5),
              "power_NW" : 0.0,
              "power_Ve" : 0.5,
              "zmin" :-12.3e-6,
              "zmax" :-10.7e-6,
              "npoints" : 1000,
              "nstates" : 20
              }
config = {
          "scan_variables" : ["power_NE"],
          "export_results" : True,
          "is_debug" : True
          }

sim = vertical_lattice_spectrum(parameters, config, output_directory=dir)
sim.runscan()
plt.show()
