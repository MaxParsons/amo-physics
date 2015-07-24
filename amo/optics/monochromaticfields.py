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
    @classmethod
    def from_lambda(cls, field_lambda):
        obj = cls()
        obj.field = field_lambda
        return obj
    
    @classmethod
    def from_combined_fields(cls, field_lambdas):
        obj = cls()
        obj.field = lambda r: np.sum(np.array([lamb(r) for lamb in field_lambdas]), axis=0)
        return obj
    
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
    def __init__(self, theta, phi, waist, power, wavelength):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis
        @param phi: Propagation angle counterclockwise from x-axis
        @param origin: location of waist

        """
        self.wavelength = wavelength
        self.wavevector = 2.0 * np.pi / wavelength
        self.waist = waist
        self.power = power
        # self.origin = origin FIXME
        self.phi = np.pi * phi / 180.0
        self.theta = np.pi * phi / 180.0
        self.e0 = 1.0
        
        
    @property
    def rayleigh_range(self):
        return np.pi * self.waist ** 2 / self.wavelength
        
    def field(self, r):
        z = r[2] * np.cos(self.theta) - r[0] * np.sin(self.phi)
        rsquared = r[1] ** 2 + (r[0] * np.cos(self.theta) + r[2] * np.sin(self.theta)) ** 2
        zr = self.rayleigh_range
        w = self.waist * np.sqrt(1 + (z / zr) ** 2)
        rad = z * (1 + (zr / z) ** 2)
        field = self.e0 * self.waist / w * np.exp(-rsquared / w ** 2 - 1.j * self.wavevector * z - 
                                                  1.j * self.wavevector * rsquared / (2. * rad) + 
                                                  1.j * np.arctan(z / zr))
        return field
        
    def get_field_lambda(self):
        lbda = lambda r: self.field(r)
        return lbda
    
if __name__ == '__main__':
    dimp = SymmetricGaussianBeam(0.0, 0.0, 50.0e-6, 0.0, 1064.0e-6)
    dimp_lambda = dimp.get_field_lambda()
    r = [0, 0, 40e-6]
    print np.abs(dimp_lambda(r)) ** 2
        
        
        
