'''
Created on Aug 21, 2015

@author: Max
'''
class TrapCalibrations(object):
    def __init__(self):
        self.NE_t_retro = 0.97
        self.NW_t_retro = 0.97
        
        self.NE_wavelength = 1.064e-6
        self.NW_wavelength = 1.064e-6
        self.Ve_wavelength = 1.064e-6
        self.accordion_wavelength = 1.064e-6
        self.dimple_wavelength = 0.78e-6
        
        self.NE_waist = 108.0e-6
        self.NW_waist = 110.1e-6
        self.Ve_waist = 42.9e-6
        self.accordion_waist = 200.0e-6
        self.dimple_waist = 30.0e-6
        
        self.NE_theta = 69.0
        self.NW_theta = 69.0
        self.Ve_theta = 5.86
        self.accordion_theta = 73.0
        self.NE_phi = 45.0
        self.NW_phi = -45.0
        self.Ve_phi = 0.0
        self.accordion_phi = 20.0
        
        self.NE_phase = 0.0
        self.NW_phase = 0.0
        
        self.NE_offset = [0.0, 0.0, 0.0]
        self.NW_offset = [0.0, 0.0, 0.0]
        self.Ve_offset = [0.0, 0.0, 0.0]
        self.dimple_offset = [0.0, 0.0, 0.0]
        self.accordion_offset = [0.0, 0.0, 0.0]
        
        self.substrate_distance = -10.995e-6# This is garbage.  I've just made sure to center on a minimum
