'''
Created on Oct 30, 2014

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
import amo.core.physicalconstants
import amo.core.misc
rot = amo.core.misc.Rotations
c = amo.core.physicalconstants.PhysicalConstantsSI

class MonochromaticField(object):
    def field(self, r, t=0):
        pass
        
    def get_amplitude(self, r, t=0, pol=None):
        if pol is None:
            return np.absolute(np.linalg.norm(self.field(self, r, t), axis=1))
        else:
            return np.absolute(np.dot(self.field(r, t), np.conjugate(pol)))
    
    def get_phase(self, r, t=0, pol=np.array([1, 0])):
        return np.angle(np.dot(self.field(r, t), np.conjugate(pol)))
    
    def get_intensity_time_averaged(self, r, pol=None):
        if pol is None:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.linalg.norm(self.field(r), axis=1)))
        else:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.dot(self.field(r), np.conjugate(pol))))

class GaussianBeam(MonochromaticField):
    def __init__(self, wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        """
        self.wavelength = wavelength
        self.k = 1.0 / wavelength
        self.waist = waist
        self.r0 = r0
        self.e0 = e0
        self.polarization = polarization / np.linalg.norm(polarization)
        self.theta = polar_angle
        self.phi = azimuthal_angle
        self.rayleigh_range = np.pi * waist ** 2 / wavelength
        
    def field(self, r, t=0):
        r = rot.Ry(rot.Rz(r, -np.pi * self.phi / 180.0), -np.pi * self.theta / 180.0)
        rsquared = np.sum(np.square(np.abs(r)), axis= -1)
        z = r[:, 2]
        w = self.waist * np.sqrt(1 + np.square(z / self.rayleigh_range))
        print rsquared
        print w
        rad_inv = np.multiply(z, np.power(np.square(z) + np.square(self.rayleigh_range), -1))
        return np.multiply(np.divide(self.e0, w), np.exp(np.divide(-rsquared, np.square(w))) - 1.j * (self.k * z + np.multiply(rsquared, np.square(rad_inv)) + np.arctan(z / self.rayleigh_range)))
    
if __name__ == '__main__':
    beam = GaussianBeam(1064e-9, 80e-6)
    xvals = np.matrix([[0, 0, 0], [0, 0, 80e-6]])
    print beam.get_intensity_time_averaged(xvals)
        
        
        
