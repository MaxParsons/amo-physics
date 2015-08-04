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
        
    def get_field_line(self, direction, center, extents, num1d=5000):
        direction = np.array(direction)
        direction = direction / np.sqrt(np.dot(direction, direction))
        center = np.array(center)
        param = np.linspace(extents[0], extents[1], num1d)
        x = center[0] + param * direction[0]
        y = center[1] + param * direction[1]
        z = center[2] + param * direction[2]
        return param, self.field(x, y, z)
    
    def plot_field_line_amplitude(self, ax, direction, center, extents, num1d=1000):
        param, field = self.get_field_line(direction, center, extents, num1d)
        ax.plot(param, np.absolute(field))
        
    def plot_field_line_intensity(self, ax, direction, center, extents, num1d=1000):
        param, field = self.get_field_line(direction, center, extents, num1d)
        ax.plot(param, 0.5 * c.epsilon0 * c.c * np.absolute(field) ** 2)
        
    def plot_field_line_phase(self, ax, direction, center, extents, num1d=1000):
        param, field = self.get_field_line(direction, center, extents, num1d)
        ax.plot(param, np.angle(field))
        
    def get_phase(self, x, y, z, t=0, pol=np.array([1, 0])):
        return np.angle(np.dot(self.field(x, y, z, t), np.conjugate(pol)))
    
    def get_intensity_time_averaged(self, x, y, z, pol=None):
        if pol is None:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.linalg.norm(self.field(x, y, z), axis=1)))
        else:
            return 0.5 * c.epsilon0 * c.c * np.square(np.absolute(np.dot(self.field(x, y, z), np.conjugate(pol))))

class SymmetricGaussianBeam(MonochromaticField):
    def __init__(self, theta, phi, waist, power, wavelength, offset=[0.0, 0.0, 0.0], phase=0.0):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis in degrees
        @param phi: Propagation angle counterclockwise from x-axis in degrees
        @param waist: beam waist
        @param power: beam power
        @param offset: location of waist
        @param phase: overall phase offset in degrees

        """
        self.wavelength = wavelength
        self.waist = waist
        self.power = power
        self.phi = np.pi * phi / 180.0
        self.theta = np.pi * theta / 180.0
        self.offset = offset
        self.phase = np.pi * phase / 180.0
        
        
    @property
    def rayleigh_range(self):
        return np.pi * self.waist ** 2 / self.wavelength
    
    @property
    def wavevector(self):
        return 2. * np.pi / self.wavelength
    
    @property
    def I0(self):
        return 2. * self.power / (np.pi * self.waist ** 2)
    
    @property
    def e0(self):
        return np.sqrt(2. * self.I0 / (c.epsilon0 * c.c))
        
    def field(self, x, y, z):
        x = x - self.offset[0]
        y = y - self.offset[1]
        z = z - self.offset[2]
        znew = x * np.sin(self.theta) * np.cos(self.phi) + y * np.sin(self.theta) * np.cos(self.phi) + z * np.cos(self.theta)
        # znew = z * np.cos(self.theta) + x * np.sin(self.theta)
        rsquared = (x * np.cos(self.theta) - z * np.cos(self.phi) * np.sin(self.theta)) ** 2 + \
            (y * np.cos(self.phi) * np.sin(self.theta) - x * np.sin(self.theta) * np.sin(self.phi)) ** 2 + \
            (y * np.cos(self.theta) - z * np.sin(self.theta) * np.sin(self.phi)) ** 2
        zr = self.rayleigh_range
        w = self.waist * np.sqrt(1 + (znew / zr) ** 2)
        rad = znew * (1 + (zr / znew) ** 2)
        field = self.e0 * self.waist / w * np.exp(-rsquared / w ** 2 - 1.j * self.wavevector * znew - 
                                                  1.j * self.wavevector * rsquared / (2. * rad) + 
                                                  1.j * np.arctan(znew / zr) - 1.j * self.phase)
        
        
        return field
    
class Lattice1d(MonochromaticField):
    def __init__(self, theta, phi, waist, power, wavelength, offset=[0.0, 0.0, 0.0], phase=0.0):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis in degrees
        @param phi: Propagation angle counterclockwise from x-axis in degrees
        @param waist: beam waist
        @param power: beam power
        @param offset: location of waist
        @param phase: overall phase offset in degrees

        """
        self.wavelength = wavelength
        self.waist = waist
        self.power = power
        self.phi = np.pi * phi / 180.0
        self.theta = np.pi * theta / 180.0
        self.offset = offset
        self.phase = np.pi * phase / 180.0
        
        self.beam0 = SymmetricGaussianBeam(theta, phi, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=self.phase)
        self.beam1 = SymmetricGaussianBeam(180.0 - theta, phi, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=self.phase + 180.0)
        
    def field(self, x, y, z):
        return self.beam0.field(x, y, z) + self.beam1.field(x, y, z)
    
class Lattice2d(MonochromaticField):
    def __init__(self, theta, phi, waist, power, wavelength, offset=[0.0, 0.0, 0.0], phase0=0.0, phase1=0.0):
        """
        __init__(wavelength, waist, polarization=np.array([1, 0]), e0=1, r0=np.array([0, 0, 0]), polar_angle=0, azimuthal_angle=0)
        all units SI, angles in degrees
        
        @param theta: Propagation angle from z-axis in degrees
        @param phi: Propagation angle counterclockwise from x-axis in degrees
        @param waist: beam waist
        @param power: beam power
        @param offset: location of waist
        @param phase: overall phase offset in degrees

        """
        self.wavelength = wavelength
        self.waist = waist
        self.power = power
        self.phi = np.pi * phi / 180.0
        self.theta = np.pi * theta / 180.0
        self.offset = offset
        self.phase0 = np.pi * phase0 / 180.0
        self.phase1 = np.pi * phase1 / 180.0
        
        self.beam0 = SymmetricGaussianBeam(theta, phi, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=phase0)
        self.beam1 = SymmetricGaussianBeam(180.0 - theta, phi, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=phase0 + 180.0)
        self.beam2 = SymmetricGaussianBeam(theta, phi + 180, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=phase1)
        self.beam3 = SymmetricGaussianBeam(180.0 - theta, phi + 180, self.waist, \
                                      self.power, self.wavelength, offset=self.offset, phase=phase1 + 180.0)
        
    def field(self, x, y, z):
        return self.beam0.field(x, y, z) + self.beam1.field(x, y, z) + self.beam2.field(x, y, z) + self.beam3.field(x, y, z)
        


def test():
    wavelength = 1.064e-6
    waist = 80.0e-6
    theta = 70.0
    phi = 0.0
    offset = [0.e-6, 0., 0.e-6]
    power = 1.0
    extents_line = [-20e-6, 20e-6]
    center_line = [0., 0., 1.e-6]
    direction = [1, 0 , 0]
    # pwave = SymmetricGaussianBeam(theta, phi, waist, 1.0, wavelength, offset=offset)
    pwave = Lattice2d(theta, phi, waist, power, wavelength, offset=offset, phase0=0.0)
    pwave_field = pwave.get_field_plane([0, 1, 0], np.array([[-20e-6, 20e-6], [0.0, 20e-6]]), 0e-6, num1d=500)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # pwave_field.plot_amplitude(ax)
    pwave.plot_field_line_phase(ax, direction, center_line, extents_line, num1d=5000)
    plt.show()
    
if __name__ == '__main__':
#    wavelength = 1.0e-6
#    pwave = PlaneWaveZTest(wavelength)
#    plt_norm = np.array([0, 0, 1])
#    plt_center = (0, 0, 0)
#    plt_extents = np.array([[-1.0e-5, 1.0e-5], [-1.0e-5, 1.0e-5]])
#    pwave.plot_intensity_plane(plt_norm, plt_center, plt_extents)
    test()
        
        
