'''
Created on Apr 6, 2015

@author: MP
'''
import numpy as np
from amo.core.physicalconstants import PhysicalConstantsSI as pc
from amo.core.physicalconstants import LithiumSixSI as li
import math
import unittest
import matplotlib.pylab as plt
from amo.quantum.montecarlowavefunction import MonteCarloWavefunction

class StochasticHeatingHarmonic(object):
    def __init__(self):
        self.frequency_gnd = 1.25*10**6
        self.frequency_ex = 1.0*10**6
        self.num_vibrational = 15
        self.rabi = 2 * np.pi * 1.0*10**5
        self.detuning = 0
        self.wavelength = pc.c / li.frequency_d1
        self.decay_rate = 2 * np.pi * 5.9 * 10**6
        self.mass = li.mass
        self.xsamples = 200
        self.tstep = 1 * 10**(-9)
        self.numsteps = 10**3
        self.sim = None
        
    @property
    def q_vibrational(self):
        return np.array(range(0,self.num_vibrational) * 2)
    
    @property
    def q_internal(self):
        return np.array[[0]*self.num_vibrational + [1] * self.num_vibrational]
    
    @property
    def k(self):
        return 2 * np.pi/self.wavelength
    
    @property
    def e_recoil(self):
        pc.hbar ** 2 * self.k**2/(2 * self.mass)
    
    def run(self):
        ham = self.get_hamiltonian()
        loss_ops = self.get_loss_operators()
        self.sim = MonteCarloWavefunction(ham, loss_ops)
        initial_state = np.zeros(2 * self.num_vibrational, dtype="complex128")
        initial_state[0] = 1.0
        self.sim.run_calculation(initial_state, self.tstep, self.numsteps)
        
    
    @property
    def xs_automatic(self):
        x0 = np.sqrt(pc.hbar/(self.mass * 2 * np.pi * np.max(np.array([self.frequency_ex, self.frequency_gnd]))))
        return np.linspace(-10 * x0, 10 * x0, self.xsamples)
    
    def get_loss_operators(self):
        length = (self.num_vibrational)*2
        nleft = np.array([range(0, self.num_vibrational) * 2,] * length)
        nright = nleft.T
        int_left = np.array([[0] * (self.num_vibrational ) + [1] * (self.num_vibrational),] * length)
        int_right = int_left.T
        compute_overlap = np.vectorize(self._compute_overlap_no_vec)
        lamb_dicke = compute_overlap(nleft, nright, int_left, int_right)
        loss_operators = []
        
        for idx_left in range(0, self.num_vibrational):
            for idx_right in range(0, self.num_vibrational):
                c = np.zeros((length, length))
                c[idx_right, self.num_vibrational + idx_left] = 1.0
                loss_operators.append(c * lamb_dicke * np.sqrt(self.decay_rate))
        
        return np.array(loss_operators)
    
    def get_hamiltonian(self):
        length = (self.num_vibrational)*2
        nleft = np.array([range(0, self.num_vibrational) * 2,] * length)
        nright = nleft.T
        int_left = np.array([[0] * (self.num_vibrational ) + [1] * (self.num_vibrational),] * length)
        int_right = int_left.T
        
        compute_h = np.vectorize(self._compute_h_el_no_vec)
        compute_overlap = np.vectorize(self._compute_overlap_no_vec)
        res = compute_h(nleft, nright, int_left, int_right) *\
         compute_overlap(nleft, nright, int_left, int_right)
        return res
          
    def _compute_h_el_no_vec(self, nleft, nright, int_left, int_right):
        if (int_left==int_right and nleft==nright):
            if int_right == 0:
                return pc.h * self.frequency_gnd * nleft
            elif int_right !=0:
                return pc.h * self.frequency_ex * nleft - pc.h * self.detuning
        elif (int_left!=int_right):
            return pc.hbar * self.rabi
        else:
            return 0
    
    def _compute_overlap_no_vec(self, nleft, nright, int_left, int_right):
        if int_left != int_right:
            xs = self.xs_automatic
            return np.abs(np.trapz(np.exp(-1.0j * self.k * xs) * self.harmonic_wavefunction(xs, nleft, int_left) *\
                 self.harmonic_wavefunction(xs, nright, int_right), xs))
        else:
            return 1.0
    
    def harmonic_wavefunction(self, x, n, int_state):
        if int_state == 0:
            f = self.frequency_gnd
        else:
            f = self.frequency_ex
        hermite_cs = np.array([0]*(n) + [1])
        x0 = np.sqrt(pc.hbar / (self.mass * 2 * np.pi * f))
        hermite = np.polynomial.hermite.hermval(x/x0, hermite_cs)
        return 1.0 / np.sqrt(2**n * math.factorial(n)) * np.sqrt(1/(np.sqrt(np.pi) * x0)) * \
            np.exp(-0.5 * (x/x0)**2) * hermite
            
    def get_nbars(self):
        return np.array([np.sum(np.conjugate(psi[0:self.num_vibrational]) * psi[0:self.num_vibrational] \
                * self.q_vibrational[0:self.num_vibrational])/ \
         np.sum(np.conjugate(psi[0:self.num_vibrational]) * psi[0:self.num_vibrational]) for psi in self.sim.psis[0:-1]])
        
    def get_motional_energies(self):
        return self.get_nbars() * pc.h * self.frequency_gnd
    
    def get_excited_population(self):
        return np.array([np.sum(np.conjugate(psi[self.num_vibrational:]) * psi[self.num_vibrational:]) for psi in self.sim.psis[0:-1]])
    
    def get_ground_populations(self):
        return np.array([np.sum(np.conjugate(psi[0:self.num_vibrational]) * psi[0:self.num_vibrational]) for psi in self.sim.psis[0:-1]])
    
    def plot_nbars(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        nbars = self.get_nbars()
        ax.plot(self.sim.ts * 10**6, nbars)
        ax.set_xlabel("time (us)")
        ax.set_ylabel("nbar")
        
    def plot_internal_populations(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pop_gnd = self.get_ground_populations()
        pop_ex = self.get_excited_population()
        ax.plot(self.sim.ts * 10**6, pop_gnd, color = "b")
        ax.plot(self.sim.ts * 10**6, pop_ex, color = "r")
        ax.set_xlabel("time (us)")
        ax.set_ylabel("internal populations")
    
        
            
class TestStochasticHeating(unittest.TestCase):
    @staticmethod
    def _test_harmonic_wavefunction():
        xs = np.linspace(-300.0*10**(-9), 300.0*10**(-9), 200)
        sh = StochasticHeatingHarmonic()
        psi = sh.harmonic_wavefunction(xs, 15, 1.25*10**6)
        print "Integral: " + str(np.trapz(psi**2, xs))
        plt.plot(xs, psi**2)
        plt.show()
    
    @staticmethod    
    def _test_compute_overlap():
        sh = StochasticHeatingHarmonic()
        ol = sh._compute_overlap_no_vec(3, 1, 0, 1)
        print "Integral: " + str(ol)
        plt.show()
        
    @staticmethod
    def _test_get_hamiltonian():
        sh = StochasticHeatingHarmonic()
        sh.num_vibrational = 3
        sh.decay_rate = 0;
        res = sh.get_hamiltonian()
        isHermitian = (np.matrix(res) - np.matrix(res).H)
        print str(res)
        print "Is Hamiltonian Hermitian?: " + str(isHermitian)
        
    @staticmethod
    def _test_loss_operators():
        sh = StochasticHeatingHarmonic()
        sh.num_vibrational = 15
        res = sh.get_loss_operators()
        
    @staticmethod
    def test_run():
        sh = StochasticHeatingHarmonic()
        sh.num_vibrational = 5
        sh.frequency_ex = sh.frequency_gnd
        sh.detuning = 0.0e6
        sh.numsteps = 20000
        sh.tstep = 2*10**(-9)
        sh.rabi = 10.0e6
        #sh.decay_rate = 0
        sh.run()
        sh.sim.plot_jumps()
        sh.plot_nbars()
        sh.plot_internal_populations()
        plt.show()
   

if __name__ == "__main__":
    unittest.main()
    plt.show()
    