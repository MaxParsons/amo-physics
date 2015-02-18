'''
Created on Feb 18, 2015

@author: Max
'''
import numpy as np
import numpy.matlib
from scipy.integrate import ode
import matplotlib.pyplot as plt

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
        

def spectrum_constant_hamiltonian():
    raman = RamanTransition()
    detunings = np.linspace(-2.0e6, 2.0e6, 100)
    four_pops = np.zeros_like(detunings)
    nbars = np.zeros_like(detunings)
    
    raman.constant_rabi = 80.0e3
    raman.anharmonicity = 26.0e3
    raman.simulation_duration = 40e-6
    raman.simulation_nsteps = 20
    raman.trap_frequency = 1.0e6
    raman.lamb_dicke = 0.28
    raman.initial_state[0] = np.sqrt(0.7)
    raman.initial_state[1] = np.sqrt(0.3)
    
    fig, ax = plt.subplots(1, 1)
    ax.set_title("simulated raman spectrum")
    ax.set_xlabel("detuning (kHz)")
    ax.set_ylabel("population in |4>")
    
    for idx, detuning in enumerate(detunings):
        print "idx = " + str(idx)
        raman.constant_detuning = detuning
        raman.compute_dynamics()
        four_pops[idx] = raman.pops_excited[-1]
        nbars[idx] = raman.nbars[-1]
    ax.plot(detunings / 1.0e3, four_pops)
    ax.plot(detunings / 1.0e3, nbars, color="k")
    plt.show()
    
        
    
if __name__ == "__main__":
    spectrum_constant_hamiltonian()
