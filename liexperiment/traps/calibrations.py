'''
Created on May 13, 2015

@author: MP
'''
import numpy as np
from plastia.core.real import RealValueErrorPair
from plastia.core.real import Real


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
    