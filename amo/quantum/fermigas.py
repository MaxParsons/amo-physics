'''
Created on Dec 11, 2014

@author: Max
'''
import amo.core.physicalconstants
import amo.core.polylog
import numpy as np
import scipy.optimize
c = amo.core.physicalconstants.PhysicalConstantsSI

class fermigasharmonic(object):
    def __init__(self, temperature, frequencies, mass, atom_number):
        self.temperature = temperature
        self.beta = 1.0 / (c.kb * temperature)
        self.frequencies = frequencies
        self.fbar = np.product(frequencies) ** (1.0 / 3.0)
        self.mass = mass
        self.atom_number
        
    def trap(self, r):
        return 0.5 * self.mass * np.dot(np.square(2.0 * np.pi * self.frequencies), np.square(r))
    
    @property
    def mu(self):
        def root_fun(mu):
            self.atom_number / (c.h * self.fbar / (c.kb * self.temperature))**3 +\
            
        scipy.optimize.root
        
