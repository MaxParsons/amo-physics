'''
Created on Aug 3, 2015

@author: Max
'''
import numpy as np
from amo.core.physicalconstants import PhysicalConstantsSI as c

def polarizability(wavelength):
    au = 1.6487772731e-41
    if wavelength == 1064.0e-9:
        return 271.0
    else:
        return None
