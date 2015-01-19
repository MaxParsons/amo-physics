'''
Created on Jan 15, 2015

@author: MP
'''
class vertical_lattice_composite(object):
    def __init__(self, f_NE, f_NW, f_VE, th_NE, th_NW, th_Ve):
        self.f_NE = f_NE
        self.f_NW = f_NW
        self.f_VE = f_VE
        self.th_NE = th_NE
        self.th_NW = th_NW
        self.th_Ve = th_Ve
        
    def depth(self, z):
        f_NW = 