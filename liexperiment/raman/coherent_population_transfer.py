'''
Created on Feb 18, 2015

@author: Max
'''
import numpy as np
import numpy.matlib
from scipy.integrate import ode
import matplotlib.pyplot as plt
from itertools import product

class RamanTransition(object):
    def __init__(self):
        self.n_vibrational = 5
        self.trap_frequency = 0.5e6
        self.anharmonicity = 26.0e3
        self.lamb_dicke = 0.28
        self.initial_state = np.zeros(2 * self.n_vibrational, dtype="complex64")
        # self.initial_state[0] = 1.0 / np.sqrt(2.0)
        self.initial_state[0] = 1.0 / np.sqrt(2.0)
        self.initial_state[1] = 1.0 / np.sqrt(2.0)
        self.constant_rabi = 500.0e3
        self.constant_detuning = -500.0e3
        self.simulation_duration = 10.0 / self.constant_rabi
        self.simulation_nsteps = 500.0
        
        nums = np.array(range(0, self.n_vibrational + 1) + range(0, self.n_vibrational + 1))
        print list(product(nums, nums))
        
        return
        
        # simulation results
        self.pops = None
        self.pops_ground = None
        self.pops_excited = None
        self.nbars = None
        self.wavefunctions = None
        self.times = None

    
    def trap_energies(self, n):
        return 2 * np.pi * (n * self.trap_frequency - 0.5 * (n - 1) * n * self.anharmonicity)
    
    def detuning(self, t):
        return 2 * np.pi * self.constant_detuning
    
    def rabi(self, t):
        return 2.0 * np.pi * self.constant_rabi
    
#    def nfactor(self, m, n):
#        if m == n:
#            return 1.0
#        elif m > n:
#            facs = np.arange(m, n)
#            return np.product(np.sqrt(facs))
#        elif m < n:
#            facs = np.arange(m, n)
#            return np.product(np.sqrt(facs + 1))
    
    def hamiltonian(self, t):
        ham0 = numpy.matlib.zeros((2 * self.n_vibrational, 2 * self.n_vibrational), dtype="complex64")
        ham1 = numpy.matlib.zeros((2 * self.n_vibrational, 2 * self.n_vibrational), dtype="complex64")
        
        for n in range(0, self.n_vibrational):
            ham0[n, n] = self.trap_energies(n)
        for n in range(self.n_vibrational, 2 * self.n_vibrational):
            ham0[n, n] = self.trap_energies(n - self.n_vibrational) - self.detuning(t)
            
        for m in range(0, self.n_vibrational):
            for n in range(self.n_vibrational, 2 * self.n_vibrational):
                ham1[m, n] = 0.5 * self.lamb_dicke ** np.abs((n - self.n_vibrational) - m) * self.rabi(t) * \
                np.exp(-1.0j * (self.detuning(t) - (self.trap_energies(n - self.n_vibrational) - self.trap_energies(m))) * t)
                       
        ham1 += ham1.H
        return ham0 + ham1
    
#    def hamiltonian(self, t):
#        ham0 = numpy.matlib.zeros((2 * self.n_vibrational, 2 * self.n_vibrational), dtype="complex64")
#        ham1 = numpy.matlib.zeros((2 * self.n_vibrational, 2 * self.n_vibrational), dtype="complex64")
        
        
    def _rhs(self, t, y):
        return 1.0j * np.dot(self.hamiltonian(t), y)
        
    def compute_dynamics(self):
        r = ode(self._rhs).set_integrator('zvode')
        r.set_initial_value(self.initial_state, 0.0)
        t1 = self.simulation_duration
        dt = t1 / self.simulation_nsteps
        ts = []
        ts.append(0.0)
        ys = []
        ys.append(self.initial_state)
        while r.successful() and r.t < t1:
            r.integrate(r.t + dt)
            ts.append(r.t)
            ys.append(r.y)
        self.times = np.array(ts)
        self.wavefunctions = np.array(ys)
        self.pops = np.abs(ys) ** 2
        self.pops_ground = np.sum(self.pops[:, 0:self.n_vibrational - 1], axis=1)
        self.pops_excited = np.sum(self.pops[:, self.n_vibrational:-1], axis=1)
        vib_states = np.append(np.arange(0, self.n_vibrational), np.arange(0, self.n_vibrational))
        self.nbars = np.sum(self.pops * vib_states, axis=1)
