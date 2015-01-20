'''
Created on Jan 19, 2015

@author: Max
'''
from datetime import datetime


class simulation(object):
    def __init__(self, variables, output_directory=None):
        self.variables = variables
        self.scan_variables = variables["scan_variables"]
        self.output_directory = output_directory
    
    @staticmethod
    def talk(instr):
        now = datetime.now()
        print '[TALK] ' + now.strftime("%H:%M:%S.%f").rstrip('0') + '--- ' + instr
        
    def print_results(self):
        pass
    
    def plot1d(self):
        pass
    
    def plot2d(self):
        pass
