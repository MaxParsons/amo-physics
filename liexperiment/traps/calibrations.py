'''
<<<<<<< HEAD
Created on May 13, 2015

@author: MP
'''
import numpy as np
from plastia.core.real import RealValueErrorPair
from plastia.core.real import Real

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

class TrapCalibrationsDipole(object):
    @staticmethod
    def crossed_dipole_trap_frequency(P_red_sheet, P_dimple, P_odt=0.0):
        f_red_x = 2.26e3 * np.sqrt(P_red_sheet)
        f_red_z = 6.78e3 * np.sqrt(P_red_sheet)
        f_dimp_r = 1.74e3 * np.sqrt(P_dimple)
        f_odt_r = 0.791e3 * np.sqrt(P_odt)
        
        f_x = np.sqrt(f_red_x**2 + f_dimp_r**2)
        f_y = np.sqrt(f_dimp_r**2 + f_odt_r**2)
        f_z = np.sqrt(f_red_z**2 + f_odt_r*2)
        f_x_err = 0
        f_y_err = 0
        f_z_err = 0
        
        fx = RealValueErrorPair(Real(f_x,unit="Hz"),Real(f_x_err, unit="Hz"))
        fy = RealValueErrorPair(Real(f_y,unit="Hz"),Real(f_y_err, unit="Hz"))
        fz = RealValueErrorPair(Real(f_z,unit="Hz"),Real(f_z_err, unit="Hz"))
        
        return fx, fy, fz
    
=======
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
>>>>>>> cd411f545b679ccac7edd8a79f3375eb1afc2042
