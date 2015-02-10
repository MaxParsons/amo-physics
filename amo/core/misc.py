'''
Created on Dec 1, 2014

@author: Max
'''
import numpy as np
class Rotations(object):
    @staticmethod
    def Rx(theta):
        """
        Rx(angle)
        
        Rotation matrix for rotating by angle (radians) about x-axis
        """
        rot = np.matrix([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        return rot
    
    @staticmethod
    def Ry(theta):
        """
        Rx(angle)
        
        Rotation matrix for rotating by angle (radians) about y-axis
        """
        rot = np.matrix([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
        return rot
    
    @staticmethod
    def Rz(theta):
        """
        Rx(angle)
        
        Rotation matrix for rotating by angle (radians) about z-axis
        """
        rot = np.matrix([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        return rot
