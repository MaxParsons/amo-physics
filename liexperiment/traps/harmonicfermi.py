'''
Created on May 22, 2015

@author: MP
'''
import amo.core.polylog as polylog
import numpy as np
import amo.core.physicalconstants
from amo.core.physicalconstants import LithiumSixSI as li
import scipy.optimize
import unittest
import mpmath
from liexperiment.traps.calibrations import TrapCalibrationsDipole as cal
pc = amo.core.physicalconstants.PhysicalConstantsSI

class HarmonicFermi(object):
    def __init__(self, frequencies, atomnumber, temperature):
        """
        @param frequencies: list of trap frequencies in Hz
        @param temperature: temperature in kelvin
        """
        self.frequencies = np.array(frequencies)
        self.atomnumber = atomnumber
        self.temperature = temperature
        self.m = li.mass
        
    @property
    def omegabar(self):
        return np.power(np.product(2 * np.pi * self.frequencies), 1.0 / self.frequencies.shape[0])
    
    @property
    def beta(self):
        return 1.0 / (pc.kb * self.temperature)
    
    @property
    def fermi_energy(self):
        return pc.hbar * (6 * self.atomnumber)**(1.0/3.0) * self.omegabar
    
    @property
    def debroglie(self):
        return np.sqrt((2 * np.pi * pc.hbar**2)/(self.m * pc.kb * self.temperature))
    
    @property
    def chemical_potential(self):
        sol = scipy.optimize.root(self._chemical_potential_eqn, np.array([self.fermi_energy]))
        return sol.x[0]

    def _chemical_potential_eqn(self, mu):
        return -self.atomnumber + polylog.fermi_poly3(self.beta * mu)*(pc.kb * self.temperature / (pc.hbar * self.omegabar))**3
    
    def central_density(self):
        return mpmath.fp.re(-1.0 / self.debroglie**3 * mpmath.polylog(3.0/2.0, -np.exp(self.beta * (self.chemical_potential))))
    
    def average_velocity(self):
        return np.sqrt(3 * pc.kb * self.temperature/self.m)
    
def compute_crossed_dipole_parameters():
    redsheetP = 0.1;
    dimpleP = 0;
    odtP = redsheetP * 11.0;
    atomnumber = 5.0e4;
    temperature = 1.0e-6;
    a = 300 * 56.0e-12;
    
    fx, fy, fz = cal.crossed_dipole_trap_frequency(redsheetP, dimpleP, odtP)
    frequencies = np.array([fx.standard_value, fy.standard_value, fz.standard_value])
    trap = HarmonicFermi(frequencies, atomnumber, temperature)
    fermir = np.sqrt(2 * trap.fermi_energy/(li.mass * (2 * np.pi * frequencies)**2))*10**6
    print "fx: {:.2E}".format(frequencies[0])
    print "fy: {:.2E}".format(frequencies[1])
    print "fz: {:.2E}".format(frequencies[2])
    print "temperature {:.2E} kelvin".format(temperature)
    print "fermi temperature: {:.2E} kelvin".format(trap.fermi_energy / pc.kb)
    print "fermi radii ({}, {}, {}) micron".format(fermir[0], fermir[1], fermir[2])
    print "Central Density: {:.2E} 1/cm^3".format(trap.central_density()/1.0e6)
    print "Central Collision Rate: {:.2E} 1/s".format(trap.average_velocity() * 4 * np.pi * a**2 * trap.central_density()) 
    
class TestHarmonicFermi(unittest.TestCase):
    def test_mu(self):
        frequencies = np.array([2.75e3, 2.75e3, 12.2])
        atomnumber = 50e4
        temperature = 5e-6
        a = 300.0
        trap = HarmonicFermi(frequencies, atomnumber, temperature)
        fermir = np.sqrt(2 * trap.fermi_energy/(li.mass * (2 * np.pi * frequencies)**2))*10**6
        print "fermi temperature: {:.2E} kelvin".format(trap.fermi_energy / pc.kb)
        print "fermi radii ({}, {}, {}) micron".format(fermir[0], fermir[1], fermir[2])
        print "Central Density: {:.2E} 1/cm^3".format(trap.central_density()/1.0e6)
        print "Central Collision Rate: {:.2E} 1/s".format(trap.average_velocity() * 4 * np.pi * a**2 * trap.central_density()) 

if __name__ == "__main__":
    unittest.main()
    #compute_crossed_dipole_parameters()
        
        