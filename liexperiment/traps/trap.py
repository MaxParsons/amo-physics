'''
Created on Aug 3, 2015

@author: Max
'''
import numpy as np
from liexperiment.lithium.polarizability import polarizability
from amo.core.physicalconstants import PhysicalConstantsSI as c
import amo.optics.monochromaticfields as field
from liexperiment.traps.calibrations import TrapCalibrations
import matplotlib.pylab as plt
import unittest
from astropy.coordinates.tests import accuracy

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
    def __init__(self, theta, phi, waist, power, wavelength, offset, phase0=0.0, visibility=1.0):
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=phase0, visibility=visibility)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        self.potential = lat.potential

class Lattice1D(Trap):
    def __init__(self, theta, phi, waist, power, wavelength, offset):
        lat_field = field.Lattice1d(theta, phi, waist, power, wavelength, offset=offset)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        self.potential = lat.potential
    
class LiLattice(Trap):
    def __init__(self, cal, NE_power, NW_power, Ve_power, accordion_power=0.0, dimple_power=0.0):
        self.NE_power = NE_power
        self.NW_power = NW_power
        self.Ve_power = Ve_power
        self.accordion_power = accordion_power
        self.dimple_power = dimple_power
        self.calibration = cal
        self.NE_spacing = cal.NE_theta
        
        traps = []
        # NE
        NE_field = field.Lattice2d(cal.NE_theta, cal.NE_phi, cal.NE_waist, NE_power, cal.NE_wavelength, offset=cal.NE_offset, \
                                   phase=cal.NE_phase, t_retro=cal.NE_t_retro)
        self.NE = Trap.from_laser_field(NE_field, polarizability(NE_field.wavelength))
        traps.append(self.NE)
        
        # NW
        NW_field = field.Lattice2d(cal.NW_theta, cal.NW_phi, cal.NW_waist, NW_power, cal.NW_wavelength, offset=cal.NW_offset, \
                                   phase=cal.NW_phase, t_retro=cal.NW_t_retro)
        self.NW = Trap.from_laser_field(NW_field, polarizability(NW_field.wavelength))
        traps.append(self.NW)
        
        # Ve
        Ve_field = field.Lattice1d(cal.Ve_theta, cal.Ve_phi, cal.Ve_waist, self.Ve_power, cal.Ve_wavelength, offset=cal.Ve_offset)
        self.Ve = Trap.from_laser_field(Ve_field, polarizability(Ve_field.wavelength))
        traps.append(self.Ve)
        
        # accordion
        accordion_field = field.Lattice1d(cal.accordion_theta, cal.accordion_phi, cal.accordion_waist, self.accordion_power, cal.accordion_wavelength, offset=cal.accordion_offset)
        self.accordion = Trap.from_laser_field(accordion_field, polarizability(accordion_field.wavelength))
        traps.append(self.accordion)
            
        if dimple_power != 0.0:
            raise Exception("Dimple not implemented yet")
        
        self.potential = lambda x, y, z: np.sum(np.array([trap.potential(x, y, z) for trap in traps]), axis=0)
    
    def chemical_potential(self, l, m, n):
        pass
        
                                        
class TestTrap(unittest.TestCase):
    def _testPolarizability(self):
        wavelength = 1.064e-6
        waist = 80.0e-6
        theta = 70.0
        phi = 0.0
        offset = [0.e-6, 0.e-6, -0.57e-6]
        power = 0.5
        extents_line = [-5e-6, 5e-6]
        center_line = [0., 0., 0.e-6]
        direction = [0, 0 , 1]
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=0.0)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        lat.plot_potential_along_line(ax, direction, center_line, extents_line, yscale=1.0 / c.h / 1.0e6, num1d=5000)
        plt.show()
        
    def _testMultiTrap(self):
        wavelength = 1.064e-6
        waist = 108.0e-6
        waist1 = 42.9e-6
        theta = 69.0
        theta_vert = 5.86
        phi = 0.0
        offset = [0.e-6, 0.e-6, 0.0e-6]
        power = 1.0
        power_vert = 0.0
        extents_line = [-20e-6, 20e-6]
        center_line = [0., 0., 0.0]
        direction = [0, 0, 1]
        t_retro = 0.97
        lat_field = field.Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=0.0, phase1=0.0, t_retro=t_retro)
        lat_field1 = field.Lattice1d(theta_vert, phi, waist1, power_vert, wavelength, offset=offset)
        lat = Trap.from_laser_field(lat_field, polarizability(lat_field.wavelength))
        lat1 = Trap.from_laser_field(lat_field1, polarizability(lat_field.wavelength))
        mytrap = Trap.from_multiple_traps([lat, lat1])
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        mytrap.plot_potential_along_line(ax, direction, center_line, extents_line, yscale=1.0 / c.h / 1.0e6, num1d=5000)
        plt.show()
        
    def testLiLattice(self):
        NE_power = 0.5
        NW_power = 0.5
        Ve_power = 0.2
        accordion_power = 0.0
        cal = TrapCalibrations()
        lat = LiLattice(cal, NE_power, NW_power, Ve_power, accordion_power=accordion_power)
        center_line = [0., 0., cal.substrate_distance]
        extents_line = [-50.e-6, 50.e-6]
        direction = [0, 1, 0]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        lat.plot_potential_along_line(ax, direction, center_line, extents_line, yscale=1.0 / c.h / 1.0e6, num1d=5000)
        plt.show()
        
if __name__ == "__main__":
    unittest.main()
    plt.show()
