'''
Created on Dec 12, 2014

@author: MP
'''
import amo.core.physicalconstants
import numpy as np
c = amo.core.physicalconstants.PhysicalConstantsSI

class latticecalibrations(object):
    recoil_NE = 25.75e3*c.h #calculated assuming 20 degree angle
    recoil_NW = 25.75e3*c.h
    recoil_Ve = 28.86e3*c.h #calculated form measured angle, lab book 8, pg. 133
    
    modulationdepth_per_watt_NE = 29.0 * recoil_NE
    modulationdepth_per_watt_NW = 30.5 * recoil_NW
    modulationdepth_per_watt_Ve = 5.2 * recoil_Ve
    
    waist_NE = 80.0e-6 #not calibrated
    waist_NW = 80.0e-6 
    waist_Ve = 80.0e-6
    
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
        depth_Er = self.depth_NE / self.cal.recoil_NE
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er *self.cal.recoil_NE / c.h
        return f
    
    @property
    def f_NW(self):
        depth_Er = self.depth_NW / self.cal.recoil_NW
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er *self.cal.recoil_NW / c.h
        return f
    
    @property
    def f_Ve(self):
        depth_Er = self.depth_Ve / self.cal.recoil_Ve
        f_Er = np.sqrt(4.0 * depth_Er)
        f = f_Er *self.cal.recoil_Ve / c.h
        return f
    
    