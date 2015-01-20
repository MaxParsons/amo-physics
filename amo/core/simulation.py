'''
Created on Jan 20, 2015

@author: Max
'''
import amo.core.utilities as util
import numpy as np
import collections
import os.path
import matplotlib.pyplot as plt

class simulation(object):
    def __init__(self, parameters, config, output_directory=None):
        self.parameters = parameters
        self.config = config
        self.scan_variables = self.config.get("scan_variables")
        self.output_directory = output_directory
        self.results = [{}]
        self.shot_figures = []
        
    def shot(self, shotnumber, **params):
        pass
    
    def runscan(self):
        if self.scan_variables is None:
            self.shot(0, self.parameters)
        else:
            scan_vars = util.cartesian([self.parameters[key] for key in self.scan_variables])
            shot_params = self.parameters.copy()
            for i, v in enumerate(scan_vars):
                for key_idx, key in enumerate(self.scan_variables):
                    shot_params[key] = v[key_idx]
                self.results.append({})
                self.results[i].update(shot_params)
                self.shot(i, shot_params)
                
        if self.config.get("export_results"):
            self.export_results()
                
                
    def export_results(self):
        if self.output_directory is None:
            raise Exception("No output directory is specified.")
        else:
            results_dir = os.path.join(self.output_directory, "results.txt")
            results_file = open(results_dir, "w")
        for key, val in self.results[0].iteritems():
            if not isinstance(val, collections.Iterable):
                results_file.write(key)
                results_file.write("\t")
        results_file.write("\n")
        for shot_result in self.results:
            for key, val in shot_result.iteritems():
                if not isinstance(val, collections.Iterable):
                    results_file.write(str(val))
                    results_file.write("\t")
            results_file.write("\n")
        results_file.close()
        
    def export_plots(self):
        pass
    
    def plot1d(self):
        pass
    
    def plot2d(self):
        pass
    
    def _generate_shot_subplots(self, ncolumns=4.0):
        nshots = 1
        for var in self.scan_variables:
            nshots = nshots * len(self.parameters[var])
            
        nrows = np.ceil(nshots / ncolumns)
        fig, axes = plt.subplots(nrows, int(ncolumns))
        axes.reshape(nrows * ncolumns,)
        return (fig, axes)
        
        
