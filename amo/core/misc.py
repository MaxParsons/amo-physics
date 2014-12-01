'''
Created on Dec 1, 2014

@author: Max
'''
import numpy as np
class Rotations(object):
    @staticmethod
    def Rx(r, theta):
        """
        Rx(vector, angle)
        
        Rotate vector (or each row of a matrix) r by theta (radians) about x-axis
        """
        rot = np.matrix([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        return np.transpose(np.dot(rot, np.transpose(r)))
    
    @staticmethod
    def Ry(r, theta):
        """
        Rx(vector, angle)
        
        Rotate vector (or each row of a matrix) r by theta (radians) about y-axis
        """
        rot = np.matrix([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
        return np.transpose(np.dot(rot, np.transpose(r)))
    
    @staticmethod
    def Rz(r, theta):
        """
        Rx(vector, angle)
        
        Rotate vector (or each row of a matrix) r by theta (radians) about z-axis
        """
        rot = np.matrix([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        return np.transpose(np.dot(rot, np.transpose(r)))
