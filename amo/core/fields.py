'''
Created on Jul 30, 2015

@author: Max
'''
import numpy as np
from colorsys import hls_to_rgb
import matplotlib.pyplot as plt

class ComplexScalarField2D(object):
    def __init__(self, coords, values):
        self.coords = coords
        self.values = values
        
    @property
    def phase(self):
        return np.angle(self.values)
    
    @property
    def amplitude(self):
        return np.absolute(self.values)
    
    @staticmethod
    def colorize(z):
        n, m = z.shape
        c = np.zeros((n, m, 3))
        c[np.isinf(z)] = (1.0, 1.0, 1.0)
        c[np.isnan(z)] = (0.5, 0.5, 0.5)
        
        idx = ~(np.isinf(z) + np.isnan(z))
        A = (np.angle(z[idx]) + np.pi) / (2 * np.pi)
        A = (A + 0.5) % 1.0
        B = 1.0 - 1.0 / (1.0 + abs(z[idx]) ** 0.3)
        c[idx] = [hls_to_rgb(a, b, 0.8) for a, b in zip(A, B)]
        return c    
    
    def plot_color_phase(self, ax):
        ax.imshow(self.colorize(self.values), interpolation='none')
        
    def plot_amplitude(self, ax):
        ax.imshow(self.amplitude)
        
