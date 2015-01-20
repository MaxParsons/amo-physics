'''
Created on Dec 12, 2014

@author: MP
'''
import amo.core.physicalconstants
import numpy as np
c = amo.core.physicalconstants.PhysicalConstantsSI
li = amo.core.physicalconstants.LithiumSixSI

class latticecalibrations(object):
    def __init__(self):
        self.wavelength = 1.064e-6
        
        # self.recoil_ax_NE = 25.75e3 * c.h# calculated assuming 20 degree angle
        # self.recoil_ax_NW = 25.75e3 * c.h
        # self.recoil_Ve = 28.86e3 * c.h# calculated form measured angle, lab book 8, pg. 133
        self.frequency_per_root_watt_rad_NE = 400.0e3
        self.frequency_per_root_watt_rad_NW = 400.0e3
        self.frequency_per_root_watt_ax_Ve = 481.27e3
        
        self.theta_NE = 20.0 * np.pi / 180.0
        self.theta_NW = 20.0 * np.pi / 180.0
        self.theta_Ve = 84.14 * np.pi / 180.0
        
        self.waist_NE = 80.0e-6# not calibrated
        self.waist_NW = 80.0e-6# not calibrated
        self.waist_Ve = 40.0e-6# not calibrated
    
    # NE lattice
    @property
    def lattice_spacing_ax_NE(self):
        return 0.5 * self.wavelength / np.sin(self.theta_NE)
    
    @property
    def lattice_spacing_ax_NW(self):
        return 0.5 * self.wavelength / np.sin(self.theta_NW)
    
    @property
    def lattice_spacing_ax_Ve(self):
        return 0.5 * self.wavelength / np.sin(self.theta_Ve)
    
    @property
    def lattice_spacing_rad_NE(self):
        return 0.5 * self.wavelength / np.cos(self.theta_NE)
    
    @property
    def lattice_spacing_rad_NW(self):
        return 0.5 * self.wavelength / np.cos(self.theta_NW)
    
    @property
    def recoil_rad_NE(self):
        return c.h ** 2 * (2.0 * self.lattice_spacing_rad_NE) ** (-2.0) / (2.0 * li.mass)
    
    @property
    def recoil_rad_NW(self):
        return c.h ** 2 * (2.0 * self.lattice_spacing_rad_NW) ** (-2.0) / (2.0 * li.mass)
    
    @property
    def recoil_ax_NE(self):
        return c.h ** 2 * (2.0 * self.lattice_spacing_ax_NE) ** (-2.0) / (2.0 * li.mass)
    
    @property
    def recoil_ax_NW(self):
        return c.h ** 2 * (2.0 * self.lattice_spacing_ax_NW) ** (-2.0) / (2.0 * li.mass)
    
    @property
    def recoil_ax_Ve(self):
        return c.h ** 2 * (2.0 * self.lattice_spacing_ax_Ve) ** (-2.0) / (2.0 * li.mass)
    
    @property
    def modulation_depth_per_watt_NE(self):
        return 0.25 * ((c.h * self.frequency_per_root_watt_rad_NE) ** 2 / self.recoil_rad_NE)
    
    @property
    def modulation_depth_per_watt_NW(self):
        return 0.25 * ((c.h * self.frequency_per_root_watt_rad_NW) ** 2 / self.recoil_rad_NW)
    
    @property
    def modulation_depth_per_watt_Ve(self):
        return 0.25 * ((c.h * self.frequency_per_root_watt_ax_Ve) ** 2 / self.recoil_ax_Ve)
        
        
    
class composite_vertical_lattice(object):
    def __init__(self):
        self.cal = latticecalibrations()
        self.power_NE = 1.0
        self.power_NW = 1.0
        self.power_Ve = 1.0
    
    def axial_depth(self, z):
        return -self.power_NE * self.cal.modulation_depth_per_watt_NE * np.sin(np.pi / self.cal.lattice_spacing_ax_NE * z) ** 2 \
            - self.power_NW * self.cal.modulation_depth_per_watt_NW * np.sin(np.pi / self.cal.lattice_spacing_ax_NW * z) ** 2 \
            - self.power_NE * self.cal.modulation_depth_per_watt_Ve * np.sin(np.pi / self.cal.lattice_spacing_ax_Ve * z) ** 2
            
    
class lattice_modulation_depth(object):
    def __init__(self, power_NE, power_NW, power_Ve):
        self.power_NE = power_NE
        self.power_NW = power_NW
        self.power_Ve = power_Ve
        self.cal = latticecalibrations
    
    @property
    def depth_NE(self):
        return self.power_NE * self.cal.modulationdepth_per_watt_NE
    
    @property
    def depth_NW(self):
        return self.power_NW * self.cal.modulationdepth_per_watt_NW
    
    @property
    def depth_Ve(self):
        return self.power_Ve * self.cal.modulationdepth_per_watt_Ve
    
    @property
    def f_NE(self):
        depth_Er = self.depth_NE / self.cal.recoil_ax_NE
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er * self.cal.recoil_ax_NE / c.h
        return f
    
    @property
    def f_NW(self):
        depth_Er = self.depth_NW / self.cal.recoil_ax_NW
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er * self.cal.recoil_ax_NW / c.h
        return f
    
    @property
    def f_Ve(self):
        depth_Er = self.depth_Ve / self.cal.recoil_Ve
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er * self.cal.recoil_Ve / c.h
        return f
    
    
