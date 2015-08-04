'''
Created on Aug 3, 2015

@author: Max
'''
import numpy as np
from liexperiment.lithium.polarizability import polarizability

class Trap(object):
    @classmethod
    def from_laser_field(cls, monochromatic_field):
        obj = cls()
        polarizability = polarizability(monochromatic_field.wavelength)
        obj.potential = lambda x, y, z: polarizability * monochromatic_field.get_intensity_time_averaged(x, y, z)
        return obj
    
class Lattice(Trap):
    def __init__(self, NW_power, NE_power, Ve_power, Dimple_power, Vertical_power):
        pass
