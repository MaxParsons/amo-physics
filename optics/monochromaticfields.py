'''
Created on Oct 30, 2014

@author: Max
'''
import core.physicalconstants as pc

class MonochromaticField(object):
    def __init__(self):
        self.wavelength = 0;
        
    def get_amplitude(self, x, y, z, t):
        pass

class GaussianBeam(MonochromaticField):