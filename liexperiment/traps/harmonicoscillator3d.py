'''
Created on Nov 27, 2014

@author: MP
'''
import numpy as np
import matplotlib.pyplot as plt
import amo.core.physicalconstants
import scipy.optimize
c = amo.core.physicalconstants.PhysicalConstantsSI

class harmonicoscillator3d(object):
    def __init__(self, frequencies, cutoffs, labels=['x', 'y', 'z']):
        self.frequencies = frequencies
        self.cutoffs = cutoffs
        self.labels = labels
        self.zero_point_energies = 0.5 * c.h * frequencies
        
    def partition_function(self, temperature):
        beta = (c.kb * temperature) ** (-1)
        return np.product(np.exp(-beta * c.h * self.frequencies / 2.0) * \
                   (1 - np.exp(-beta * c.h * self.frequencies)) ** (-1))
        
    def average_energy(self, temperature):
        beta = (c.kb * temperature) ** (-1)
        return np.sum(c.h * self.frequencies * \
               (np.exp(beta * c.h * self.frequencies) - 1) ** (-1))
        
    def temperature(self, energy):
        def _average_energy_zero(temperature):
            return self.average_energy(temperature) - energy
        return scipy.optimize.brentq(_average_energy_zero, 0.1 * energy/(c.kb), 2.0 * energy/(c.kb))
        
    def population(self, states, temperature):
        beta = (c.kb * temperature) ** (-1)
        return np.product(np.exp(-beta * c.h * self.frequencies * (states + 0.5)))\
            / self.partition_function(temperature)
    
    def population_sum_over_first_frequency(self, states, temperature):
        return np.sum([self.population(np.array([n, states[0], states[1]]), temperature) for n in range(0, self.cutoffs[0])])
        pass
    
    def population_sum_over_all_but_first_frequency(self, state, temperature):
        sum = 0.0
        for n1 in np.arange(0, self.cutoffs[1]):
            for n2 in np.arange(0, self.cutoffs[2]):
                sum += self.population(np.array([state, n1, n2]), temperature)
        return sum
    
    def populations_sum_over_all_but_first_frequency(self, temperature):
        pops = np.zeros((self.cutoffs[0],))
        for n in np.arange(0, self.cutoffs[0]):
            pops[n] = self.population_sum_over_all_but_first_frequency(n, temperature)
        return pops
    
    def atoms_remaining(self, temperature):
        atoms = 0
        for idx_z in range(0, self.cutoffs[2]):
            for idx_y in range(0, self.cutoffs[1]):
                atoms += self.population_sum_over_first_frequency(np.array([idx_y, idx_z]), temperature)
        return atoms
                
                
if __name__ == '__main__':
    Er = 75.0e3 * c.h
    Tr = Er / c.kb
    frequencies = np.array([650.0e3, 950.0e3, 1300.0e3])
    cutoffs = np.array([49, 7, 11])   
    lat = harmonicoscillator3d(frequencies, cutoffs)
    temperatures = np.linspace(0.01 * Tr, 30 * Tr)
    energies = [lat.average_energy(temp) for temp in temperatures]
    populations = [lat.population(np.array([0, 0, 0]), temp) for temp in temperatures]
    radial_pops = [lat.population_sum_over_first_frequency(np.array([0, 0]), temp) for temp in temperatures]
    atoms_remaining = [lat.atoms_remaining(temp) for temp in temperatures]
    energies = np.array(energies)
     
    #plt.plot(energies / Er, atoms_remaining, marker='o')
    plt.plot(energies / Er, c.kb * temperatures / Er)
    print lat.temperature(20.0 * Er)*c.kb / Er
    plt.show()
        
