'''
Created on Jun 4, 2014

@author: Max
'''
import numpy as np
import amo.core.physicalconstants
import amo.core.utilities
misc = amo.core.utilities.simulations
import scipy.linalg
from scipy.optimize import brentq
import matplotlib.pylab as plt

pc = amo.core.physicalconstants.PhysicalConstantsSI

class Trap1D(object):
    def __init__(self, x_values, potential, mass, n_states, statistics='fermi', degeneracy=0.1, npol=200):
        self.mass = mass
        self.n_states = n_states
        self.potential = potential
        self.statistics = statistics
        self.degeneracy = degeneracy
        self.npol = npol
        self.x_values = x_values
        self.compute_eigenstates()
        
        
    def compute_eigenstates(self):# FIXME: automatic computation of cutoffs for eigenstate calculation
        misc.talk('Diaganolizing trap hamiltonian...')
        h0 = self.potential + 2.0 * (pc.hbar ** 2 / (2.0 * self.mass)) * np.ones_like(self.x_values) / (self.x_values[1] - self.x_values[0]) ** 2
        h1 = -(pc.hbar ** 2 / (2.0 * self.mass)) * np.ones_like(self.x_values) / (self.x_values[1] - self.x_values[0]) ** 2
        h_band = np.array([h0, h1])
        self.energies, self.eigenStates = scipy.linalg.eig_banded(h_band, lower=True)
        print self.energies
        idx = self.energies.argsort()
        self.energies = self.energies[idx]
        self.eigenStates = self.eigenStates[:, idx]
        print self.energies
        self.energies = self.energies[0:self.n_states]
        self.eigenStates = self.eigenStates[:, 0:self.n_states]
        self.normalize_states()
        misc.talk('Done!')
        
    def compute_chemical_potential(self):
        self.temperature = self.degeneracy * self.energies[self.npol] / pc.kb
        def chem_pot_root(mu):
            if self.statistics is 'fermi':
                return np.sum(1.0 / (np.exp((self.energies - mu) / (pc.kb * self.temperature)) + 1)) - self.npol
            elif self.statistics is 'bose':
                return np.sum(1.0 / (np.exp((self.energies - mu) / (pc.kb * self.temperature)) - 1)) - self.npol
            else:
                raise NotImplementedError('Statistics have to be "bose" or "fermi"')
                
        tolerance = self.energies[self.npol + 1] - self.energies[self.npol]
        self.chemical_potential = brentq(chem_pot_root, np.min(self.energies), np.max(self.energies), xtol=tolerance)# Use root-finding to compute chemical potential
    
    def compute_occupations(self):
        misc.talk('Computing occupations')
        self.compute_chemical_potential()
        if self.statistics is 'fermi':
            self.occupations = 1.0 / (np.exp((self.energies - self.chemical_potential) / (pc.kb * self.temperature)) + 1)
        elif self.statistics is 'bose':
            self.occupations = 1.0 / (np.exp((self.energies - self.chemical_potential) / (pc.kb * self.temperature)) - 1)
        else: 
            raise NotImplementedError('Statistics have to be bose or fermi')
        self.occupations = self.occupations / np.sum(self.occupations)
        
    def compute_density(self):
        misc.talk('Computing density distribution')
        self.compute_occupations()
        self.density = np.zeros_like(self.eigenStates[:, 0])
        for i in range(0, self.n_states):
            self.density += self.occupations[i] * np.abs(self.eigenStates[:, i]) ** 2
        self.density = self.npol * self.density
        misc.talk('Done!')
        
        
    @staticmethod
    def normalize_state(vector, xvalues):
        return vector / np.sqrt(np.trapz(np.absolute(vector) ** 2, xvalues))
    
    def normalize_states(self):
        self.eigenStates = np.apply_along_axis(self.normalize_state, 0, self.eigenStates, self.x_values)
        
    def compute_overlap(self, start_state, end_state):
        overlap = np.trapz(np.conjugate(self.eigenStates[:, start_state]) * self.eigenStates[:, end_state], self.x_values)
        de = self.energies[end_state] - self.energies[start_state] 
        return overlap, de
        
    def compute_vibrational_coupling(self, wavevector, start_state):# FIXME: vibrational coupling between some state and all others
        misc.talk('Computing vibrational couplings...')
        dip_op = np.exp(1.j * wavevector * self.x_values)
        dip_op = dip_op.reshape(len(self.x_values), 1)
        energyRescale = 1 / (pc.h * 1.e3)
        mat = np.conjugate(self.eigenStates[:, start_state].reshape(len(self.x_values), 1)) * dip_op * self.eigenStates
        trans_el = np.apply_along_axis(np.trapz, 0, mat, self.x_values)
        de = self.energies - self.energies[start_state]
        plt.xlabel('Energy (kHz)')
        plt.ylabel('Raman coupling (arb)')
        plt.plot(de * energyRescale, np.abs(trans_el))
        misc.talk('Done')
    
    def simulate_spectrum(self, wavevector):
        misc.talk('Computing spectrum...')
        dip_op = np.exp(1.j * wavevector * self.x_values)
        dip_op = dip_op.reshape(len(self.x_values), 1)
        energyRescale = 1 / (pc.h * 1.e3)
        des = []
        trans_els = []
        occupations = []
        for start_state in range(0, self.n_states):
            mat = np.conjugate(self.eigenStates[:, start_state].reshape(len(self.x_values), 1)) * dip_op * self.eigenStates
            trans_el = np.apply_along_axis(np.trapz, 0, mat, self.x_values)
            de = self.energies - self.energies[start_state]
            trans_els = np.append(trans_els, trans_el)
            occupations = np.append(occupations, self.occupations[start_state] * np.ones_like(de))
            des = np.append(des, de)
            print start_state
            
        des = np.asarray(des)
        trans_els = np.asarray(trans_els)
        print trans_els
        plt.xlabel('Energy (kHz)')
        plt.ylabel('Raman spectrum')
        plt.plot(des * energyRescale, occupations * np.abs(trans_els), '.')
        misc.talk('Done')
            
    def plot_eigenstates(self, skip=False, is_overlay_potential=True, energy_rescale=None, x_rescale=None, wavefunction_scale=0.3):
        plt.xlabel('Distance (um)')
        plt.ylabel('Energy (MHz), Re(wavefunction)')
        plt.title('Trap Eigenstates')
        if energy_rescale is None:
            energy_rescale = 1.0
        if x_rescale is None:
            x_rescale = 1.0
        if not skip:
            for i in range(0, self.n_states):
                plt.plot(self.x_values * x_rescale, np.real(self.eigenStates[:, i]) * (wavefunction_scale * energy_rescale / self.n_states * (np.max(self.energies) - np.min(self.energies)) / np.max(np.abs(self.eigenStates[:, i])))
                         + self.energies[i] * energy_rescale)
        else:
            for i in range(0, self.n_states, 10):
                plt.plot(self.x_values * x_rescale, np.real(self.eigenStates[:, i]) * (wavefunction_scale * energy_rescale / self.n_states * (np.max(self.energies) - np.min(self.energies)) / np.max(np.abs(self.eigenStates[:, i])))
                         + self.energies[i] * energy_rescale)
                
        if is_overlay_potential:
            plt.plot(self.x_values * x_rescale, self.potential * energy_rescale, color='k', linewidth=2.0)
            
    
    def plot_density_of_states(self):# FIXME: add density of states
        plt.xlabel('Energy (kHz)')
        plt.ylabel('Number of States')
        energyRescale = 1 / (pc.h * 1.e3)
        plt.hist(self.energies * energyRescale, bins=np.ceil(self.n_states / 5.0), normed=1 / self.n_states)
        
    def plot_density_lattice(self, lattice_spacing=None):
        self.compute_density()
        plt.title('Trap Density')
        if lattice_spacing is None:
            xRescale = 1. / (1.e-6)
            plt.xlabel('Distance (um)')
            plt.ylabel('Density (atoms/um')
        else:
            xRescale = 1. / (lattice_spacing)
            plt.xlabel('Distance (sites)')
            plt.ylabel('Density (atoms/site)')
            
        plt.plot(self.x_values * xRescale, self.density / xRescale)
        
    def plot_occupations(self):
        self.compute_occupations()
        plt.title('occupation fraction at {0:.2F} TF'.format(self.degeneracy))
        plt.xlabel('energy (kHz)')
        plt.ylabel('occupation fraction')
        plt.plot(self.energies / pc.h / 1e3, self.occupations)
        
        
