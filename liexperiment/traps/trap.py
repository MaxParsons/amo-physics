'''
Created on Aug 3, 2015

@author: Max
'''
import numpy as np
from liexperiment.lithium.polarizability import polarizability
from amo.core.physicalconstants import PhysicalConstantsSI as c
import amo.optics.monochromaticfields as field
import matplotlib.pylab as plt
import unittest

class Trap(object):
    @classmethod
    def from_laser_field(cls, monochromatic_field, polarizability):
        obj = cls()
        obj.potential = lambda x, y, z:-polarizability / (2 * c.epsilon0 * c.c) * monochromatic_field.get_intensity_time_averaged(x, y, z)
        return obj
    
    @classmethod
    def from_multiple_traps(cls, traps):
        obj = cls()
        obj.potential = lambda x, y, z: np.sum(np.array([trap.potential(x, y, z) for trap in traps]), axis=0)
        return obj
    
    def plot_potential_along_line(self, ax, direction, center, extents, yscale=1.0, num1d=5):
        param, potential = self.get_potential_along_line(direction, center, extents, num1d)
        ax.plot(param, yscale * potential)
    
    def get_potential_along_line(self, direction, center, extents, num1d=5000):
        direction = np.array(direction)
        direction = direction / np.sqrt(np.dot(direction, direction))
        center = np.array(center)
        param = np.linspace(extents[0], extents[1], num1d)
        x = center[0] + param * direction[0]
        y = center[1] + param * direction[1]
        z = center[2] + param * direction[2]
        return param, self.potential(x, y, z)
    
class Lattice2D(Trap):
    def __init__(self, theta, phi, waist, power, wavelength, offset, phase0=0.0):
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        self.potential = lat.potential

class Lattice1D(Trap):
    def __init__(self, theta, phi, waist, power, wavelength, offset):
        lat_field = field.Lattice1d(theta, phi, waist, power, wavelength, offset=offset)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        self.potential = lat.potential
    
class LiLattice(Trap):
    def __init__(self, NE_power, NW_power, Ve_power, accordion_power=0.0, dimple_power=0.0):
        self.NE
    
class Lattice(Trap):
    def __init__(self, NW_power, NE_power, Ve_power, Dimple_power, Vertical_power):
        pass
    
class TestTrap(unittest.TestCase):
    def _testPolarizability(self):
        wavelength = 1.064e-6
        waist = 80.0e-6
        theta = 70.0
        phi = 0.0
        offset = [0.e-6, 0.e-6, -0.57e-6]
        power = 0.5
        extents_line = [-80e-6, 80e-6]
        center_line = [0., 0., 0.e-6]
        direction = [0, 0 , 1]
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=0.0)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        lat.plot_potential_along_line(ax, direction, center_line, extents_line, yscale=1.0 / c.h / 1.0e6, num1d=5000)
        plt.show()
        
    def testMultiTrap(self):
        wavelength = 1.064e-6
        waist = 80.0e-6
        waist1 = 40.0e-6
        theta = 70.0
        theta_vert = 5.6
        phi = 0.0
        offset = [0.e-6, 0.e-6, -0.57e-6]
        power = 10
        extents_line = [-20e-6, 20e-6]
        center_line = [0., 0., 0.e-6]
        direction = [0, 0 , 1]
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=0.0)
        lat_field1 = field.Lattice1d(theta_vert, phi, waist1, power, wavelength, offset=offset)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        lat1 = Trap.from_laser_field(lat_field1, polarizability(lat_field.wavelength))
        mytrap = Trap.from_multiple_traps([lat, lat1])
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        mytrap.plot_potential_along_line(ax, direction, center_line, extents_line, yscale=1.0 / c.h / 1.0e6, num1d=5000)
        plt.show()
        
if __name__ == "__main__":
    unittest.main()
    plt.show()
