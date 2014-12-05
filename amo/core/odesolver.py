'''
Created on Dec 5, 2014

@author: Max
'''
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

class rateequation(object):
    def __init__(self, transition_matrix, initial_populations):
        """
        Initialize a rate equation solver.

        @param transition_matrix: transition matrix (dp/dt = (transition_matrix).p)
        @type transition_matrix: ndarray (numpy)
        @param initial_populations: populations at t=0
        @type initial_populations: ndarray (numpy)
        """
        self.transition_matrix = transition_matrix
        self.initial_populations = initial_populations
        
    def f(self, t, population):
        return np.dot(self.transition_matrix, population)
    
    def solve(self, duration, numsteps=10):
        times = np.linspace(0, duration, numsteps)
        populations = np.zeros((times.shape[0], self.initial_populations.shape[0]))
        populations[0, :] = self.initial_populations
        self.solver = ode(self.f)
        self.solver.set_initial_value(self.initial_populations, 0)
        for i, v in enumerate(times[1:]):
            self.solver.integrate(v)
            populations[i + 1, :] = self.solver.y
            
        self.result = rateequation._result(times, populations)
            
    class _result(object):
        def __init__(self, t, populations):
            self.t = t
            self.populations = populations
            
        def plot(self, ncutoff=None):
            fig, ax = plt.subplots(1, 1)
            if ncutoff is None:
                ncutoff = self.populations.shape[1] - 1
            for i in range(0, ncutoff + 1):
                ax.plot(self.t, self.populations[:, i], marker="o", linestyle="None", label=str(i))
                
            ax.legend()
            
if __name__ == '__main__':
    trans = np.array([[-1, 0 ], [0, -1]])
    initial_pop = np.array([1, 0])
    r = rateequation(trans, initial_pop)
    r.solve(5, numsteps=50)
    r.result.plot()
    plt.show()
