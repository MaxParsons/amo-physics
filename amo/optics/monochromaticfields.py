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

class SymmetricGaussianBeam(MonochromaticField):
    def __init__(self, k, origin, waist, power, wavelength, polarization=np.array([1, 0])):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        """
        self.wavelength = wavelength
        self.knorm = k / np.linalg.norm(k)
        self.waist = waist
        self.power = power
        self.origin = origin
        self.polarization = polarization / np.linalg.norm(polarization)
        
    @property
    def rayleigh_range(self):
        np.pi * self.waist ** 2 / self.wavelength
        
    def field(self, x, y, z, t=0):
        x, y, z = self.rotated_coordinates(x, y, z)
        rsquared = np.square(x) + np.square(y)
        w = self.waist * np.sqrt(1 + np.square(z / self.rayleigh_range))
        rad_inv = np.multiply(z, np.power(np.square(z) + np.square(self.rayleigh_range), -1))
        return np.multiply(np.divide(self.e0, w), np.exp(np.divide(-rsquared, np.square(w))) - 1.j * (self.k * z + np.multiply(rsquared, np.square(rad_inv)) + np.arctan(z / self.rayleigh_range)))
    
    def rotated_coordinates(self, x, y, z):
        phi = np.arctan(self.knorm[1] / self.knorm[0])
        theta = np.arccos(self.knorm[2])
        r = np.array([x, y, z])
        print r.shape
        print r
        rrot = np.dot(np.dot(rot.Ry(-phi), rot.Rz(-theta)), r)
        
        return rrot[0], rrot[1], rrot[2]
    
if __name__ == '__main__':
    k = np.array([1, 0, 0])
    origin = np.array([0, 0, 0])
    waist = 50.0e-6
    power = 1.0e-3
    wavelength = 1.064e-6
    beam = SymmetricGaussianBeam(k, origin, waist, power, wavelength)
    xaxis = np.linspace(0, 1.0, 5)
    yaxis = np.linspace(0, 1.0, 5)
    zaxis = np.linspace(0, 1.0, 5)
    x, y, z = np.meshgrid(xaxis, yaxis, zaxis)
    print beam.rotated_coordinates(x, y, z)
        
        
        
