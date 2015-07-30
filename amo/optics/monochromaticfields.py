'''
Created on Oct 30, 2014

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
from amo.core.fields import ComplexScalarField2D
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
    
    def get_field_plane(self, normal, extents, constant=0, num1d=100):
        a0_raw = np.linspace(extents[0, 0], extents[0, 1], num1d)
        a1_raw = np.linspace(extents[1, 0], extents[1, 1], num1d)
        a0, a1 = np.meshgrid(a0_raw, a1_raw)
        if normal[0] != 0:
            x = constant
            y = a0
            z = a1
        else:
            x = a0
            if normal[1] != 0:
                y = constant
                z = a1
            else:
                y = a1
                z = constant
        
        field = ComplexScalarField2D([x, y, z], self.field(x, y, z))
        return field
        
    def get_field_line(self, direction, extents, num1d=1000):
        pass
        
    def get_field_lambda(self):
        lbda = lambda x, y, z: self.field(x, y, z)
        return lbda
        
    def get_amplitude(self, x, y, z, t=0, pol=None):
        if pol is None:
            return np.absolute(np.linalg.norm(self.field(self, x, y, z, t), axis=1))
        else:
            return np.absolute(np.dot(self.field(x, y, z, t), np.conjugate(pol)))
    
    def get_phase(self, x, y, z, t=0, pol=np.array([1, 0])):
        return np.angle(np.dot(self.field(x, y, z, t), np.conjugate(pol)))
    
    def get_intensity_time_averaged(self, x, y, z, pol=None):
        if pol is None:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.linalg.norm(self.field(x, y, z), axis=1)))
        else:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.dot(self.field(x, y, z), np.conjugate(pol))))
    
    def plot_phase_plane(self, normal, center, extents):
        pass
    
    def plot_phase_line(self, point0, point1):
        pass

        
    def plot_intensity_plane(self, normal, center, extents):
        pass
    
    def plot_intensity_line(self, point0, point1):
        pass

class SymmetricGaussianBeam(MonochromaticField):
    def __init__(self, theta, phi, waist, power, wavelength):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis in degrees
        @param phi: Propagation angle counterclockwise from x-axis in degrees
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
        
    def field(self, x, y, z):
        znew = z * np.cos(self.theta) - x * np.sin(self.theta)
        rsquared = y ** 2 + (x * np.cos(self.theta) + z * np.sin(self.theta)) ** 2
        zr = self.rayleigh_range
        w = self.waist * np.sqrt(1 + (znew / zr) ** 2)
        rad = znew * (1 + (zr / znew) ** 2)
        field = self.e0 * self.waist / w * np.exp(-rsquared / w ** 2 - 1.j * self.wavevector * znew - 
                                                  1.j * self.wavevector * rsquared / (2. * rad) + 
                                                  1.j * np.arctan(znew / zr))
        return field
    
class PlaneWaveZTest(MonochromaticField):
    def __init__(self, wavelength):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis
        @param phi: Propagation angle counterclockwise from x-axis
        @param origin: location of waist

        """
        self.wavelength = wavelength
        self.wavevector = 2.0 * np.pi / wavelength
        self.e0 = 1.0
        
        
    @property
    def rayleigh_range(self):
        return np.pi * self.waist ** 2 / self.wavelength
        
    def field(self, x, y, z):
        field = self.e0 * np.exp(1.j * self.wavevector * z)
        return field
        
    
def func(x, y):
    return x * y

def test():
    wavelength = 1.064e-6
    waist = 2.0e-6
    theta = 45.0
    phi = 0.0
    pwave = SymmetricGaussianBeam(theta, phi, waist, 1.0, wavelength)
    pwave_field = pwave.get_field_plane([0, 1, 0], np.array([[-20e-6, 20e-6], [-20e-6, 20e-6]]), 0, num1d=500)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    pwave_field.plot_color_phase(ax)
    plt.show()
    
if __name__ == '__main__':
#    wavelength = 1.0e-6
#    pwave = PlaneWaveZTest(wavelength)
#    plt_norm = np.array([0, 0, 1])
#    plt_center = (0, 0, 0)
#    plt_extents = np.array([[-1.0e-5, 1.0e-5], [-1.0e-5, 1.0e-5]])
#    pwave.plot_intensity_plane(plt_norm, plt_center, plt_extents)
    test()
        
        
