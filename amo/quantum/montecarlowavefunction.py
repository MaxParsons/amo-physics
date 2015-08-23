'''
Created on Apr 6, 2015

@author: MP
'''
import numpy as np
from amo.core.physicalconstants import PhysicalConstantsSI as pc
import matplotlib.pyplot as plt

class MonteCarloWavefunction(object):
    def __init__(self, hamiltonian, relaxation_operators):
        """
        @param hamiltonian: hermitian part of hamiltonian
        @type hamiltonian: ndarray
        @param relaxation_operators: c_{m} operators of Molmer et. al.
        @type relaxation_operators: list of ndarrays
        """
        self.hamiltonian = hamiltonian
        self.relaxation_operators = relaxation_operators
        self.ts = None
        self.psis = None
        self.jumps = None
        self.jumps_total = None
        self.jump_probabilities = None
        
    def run_calculation(self, initial_state, timestep, nsteps):
        dt = timestep
        p_ops = [np.dot(np.conjugate(np.transpose(c)), c) for c in self.relaxation_operators]
        h = self.hamiltonian - 0.5j * pc.hbar * np.sum(np.array(p_ops), axis=0)
        h = h.astype("complex128")
        self.ts = np.linspace(0, dt * nsteps, nsteps)
        self.jumps = np.zeros((self.ts.shape[0], len(p_ops)))
        self.jump_probabilities = np.zeros_like(self.ts)
        self.jumps_total = np.zeros_like(self.ts)
        self.psis = np.zeros((nsteps+1, h.shape[0]), dtype="complex128")
        self.psis[0] = initial_state
        evolution = np.identity(initial_state.shape[0], dtype="complex128") - 1.0j/pc.hbar * timestep * h
        
        #random numbers for monte carlo
        r0 = np.random.random_sample(nsteps)
        r1 = np.random.random_sample(nsteps)
        
        njumps = 0
        for step in range(0, nsteps-1):
            self.psis[step+1] = np.dot(evolution, self.psis[step])
            dps = np.array([dt * np.vdot(self.psis[step], np.dot(p_op, self.psis[step])) for p_op in p_ops])
            dp = np.real(np.sum(dps))
            self.jump_probabilities[step+1] = dp
            
            if dp < r0[step]: #no jump
                self.psis[step+1] = self.psis[step+1] / np.linalg.norm(self.psis[step+1])
            else:  # Determine based on jump probabilities which jump to take
                print "Jump!"
                njumps+=1
                probs = dps / dp
                upto = 0
                for idx, prob in enumerate(probs):
                    if upto + prob > r1[step]:
                        step_jump = self.relaxation_operators[idx]
                        step_prob = prob
                        self.jumps[step + 1, idx] = self.jumps[step, idx] + 1 #add a jump
                    else:
                        upto += prob
                        
                self.psis[step+1] = np.dot(step_jump, self.psis[step]) \
                / np.sqrt(step_prob / dt)
                
            self.jumps_total[step+1] = njumps
                
                    
                
                
    
    def get_operator_expectation(self, operator):
        pass
    
    def get_jump_rate(self):
        pass
    
    def plot_jumps(self):
        fig = plt.figure()
        ax0 = fig.add_subplot(121)
        ax0.plot(self.ts * 10**6, self.jumps_total)
        ax0.set_xlabel("time (us)")
        ax0.set_ylabel("Number of jumps")
        
        ax1 = fig.add_subplot(122)
        ax1.plot(self.ts * 10**6, self.jump_probabilities)
        ax1.set_xlabel("time (us)")
        ax1.set_ylabel("Jump probability")
        
        
    
    def plot_jump_probabilities(self):
        pass