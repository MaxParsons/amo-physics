'''
Created on Dec 5, 2014

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
from amo.core import odesolver
from amo.core import physicalconstants
from liexperiment.traps import harmonicoscillator3d
c = physicalconstants.PhysicalConstantsSI

Er = c.h * 75.0e3

class RamanSimulation(object):
    def __init__(self):
        # SI units everywhere
        # experiment parameters
        self.f_trap = 0.65e6  # trap frequency
        self.detuning = -self.f_trap
        self.cutoff_state = 8.0  # first state to tunnel, indexed from 0
        self.pi_time = 10.0e-6  # raman carrier pi-time
        self.pump_time = 10.0e-6  # pumping time
        self.photons_per_raman = 2.0  # recoil heating for fudge factor
        self.duration = 0.5e-3  # how long is the raman cooling
        self.numsteps = 500  # how many points for the simulation
    
        # constants
        self.f_recoil = 75.0e3
        self.f_anharmonicity = 26.0e3 
    
    
    # calculated from experiment parameters
    @property
    def cutoff_state(self):
        return self._cutoff_state
    @cutoff_state.setter
    def cutoff_state(self, value):
        self._cutoff_state = value
        self.init_populations = np.zeros((self._cutoff_state + 1,))
        self.transition_matrix = np.zeros((self._cutoff_state + 1, self._cutoff_state + 1))
        
        
        
    @property
    def eta_0(self):
        return np.sqrt(self.f_recoil / self.f_trap)
    @property
    def eta(self):
        return 0.19 * np.sqrt(1.0e6 / self.f_trap)
    @property
    def eps(self):
        return self.photons_per_raman * self.eta_0 ** 2
    
    def lineshape(self, detuning, raman_rabi, decay_rate):
        s0 = 2.0 * raman_rabi ** 2 / decay_rate ** 2
        return decay_rate / 2.0 * s0 / (1 + s0 + (2.0 * detuning / decay_rate) ** 2)
    
    def sideband_frequencies(self, n0, nf):
        if n0 < nf:
            return self.f_trap * (nf - n0) - np.sign(nf - n0) * self.f_anharmonicity * np.sum(np.arange(n0, nf))
        else:
            return self.f_trap * (nf - n0) - np.sign(nf - n0) * self.f_anharmonicity * np.sum(np.arange(nf, n0))
    
    def transition_rate_neighbor_only(self, n0, nrate):
        decay_rate = 2.0 * np.pi / self.pump_time
        if nrate == self._cutoff_state:
            if (nrate - n0) == 1:  # heating from state below
                up_rabi = np.sqrt(n0 + 1) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_heat = self.sideband_frequencies(n0, n0 + 1)
                return (1 + self.eps) * (self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) + self.eta ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            elif nrate == n0: 
                return 0.0
            else:
                return 0.0
        elif nrate == self._cutoff_state - 1:
            if (nrate - n0) == 1:
                up_rabi = np.sqrt(n0 + 1) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_heat = self.sideband_frequencies(n0 , n0 + 1)
                return (1 + self.eps) * (self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) + self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            elif (nrate - n0) == -1:
                return 0.0
            elif nrate == n0:
                down_rabi = np.sqrt(n0) * self.eta * np.pi / self.pi_time
                up_rabi = np.sqrt(n0 + 1) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_heat = self.sideband_frequencies(n0, n0 + 1)
                f_cool = self.sideband_frequencies(n0, n0 - 1)
                return (1 + self.eps) * (-self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) - self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate)) + \
                        (1 - self.eps) * (-self.lineshape(self.detuning - f_cool, down_rabi, decay_rate) - self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            else:
                return 0.0
        else:
            if (nrate - n0) == 1:
                up_rabi = np.sqrt(n0 + 1) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_heat = self.sideband_frequencies(n0 , n0 + 1)
                return (1 + self.eps) * (self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) + self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            elif (nrate - n0) == -1:
                down_rabi = np.sqrt(n0) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_cool = self.sideband_frequencies(n0, n0 - 1)
                return (1 - self.eps) * (self.lineshape(self.detuning - f_cool, down_rabi, decay_rate) + self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            elif nrate == n0:
                down_rabi = np.sqrt(n0) * self.eta * np.pi / self.pi_time
                up_rabi = np.sqrt(n0 + 1) * self.eta * np.pi / self.pi_time
                stay_rabi = np.pi / self.pi_time
                f_heat = self.sideband_frequencies(n0, n0 + 1)
                f_cool = self.sideband_frequencies(n0, n0 - 1)
                if n0 == 0:
                    return (1 + self.eps) * (-self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) - self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
                else:
                    return (1 + self.eps) * (-self.lineshape(self.detuning - f_heat, up_rabi, decay_rate) - self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate)) + \
                        (1 - self.eps) * (-self.lineshape(self.detuning - f_cool, down_rabi, decay_rate) - self.eta_0 ** 2 * self.lineshape(self.detuning, stay_rabi, decay_rate))
            else:
                return 0.0
    
    def set_transition_matrix(self):
        for n0 in np.arange(0, len(self.init_populations)):
            for nrate in np.arange(0, len(self.init_populations)):
                self.transition_matrix[nrate, n0] = self.transition_rate_neighbor_only(n0, nrate)
                
    def set_mean_phonon_number(self):
        numbers = np.arange(0, self._cutoff_state + 1)
        self.mean_phonon_number = np.dot(self.populations, numbers)
        
    def simulate_cooling(self, isPlot=False):
        self.set_transition_matrix()
    
        # print np.sum(self.transition_matrix, 0)
        r = odesolver.rateequation(self.transition_matrix, self.init_populations)
        r.solve(self.duration, numsteps=self.numsteps)
        self.populations = r.result.populations
        self.times = r.result.t
        self.set_mean_phonon_number()
        if isPlot:
            r.result.plot()

def set_default_simulation_parameters(sim):
    sim.f_trap = 0.65e6  # trap frequency
    sim.f_recoil = 75e3
    sim.f_anharmonicity = 3.5e3
    sim.cutoff_state = 15.0  # first state to tunnel, indexed from 0
    trap_frequencies = np.array([sim.f_trap, 1000.0e3, 1300.0e3])
    cutoffs = np.array([sim._cutoff_state, 8, 11])   
    sim.detuning = -sim.f_trap
    sim.pi_time = 2.0e-6  # raman carrier pi-time
    sim.pump_time = 10.0e-6  # pumping time
    sim.photons_per_raman = 2.0  # recoil heating for fudge factor
    sim.duration = 0.2e-3  # how long is the raman cooling
    sim.numsteps = 100  # how many points for the simulation
    
        
if __name__ == "__main__":
    

    
    figmean, axmean = plt.subplots(1, 1)
    figgnd, axgnd = plt.subplots(1, 1)
    figlost, axlost = plt.subplots(1, 1)
    
    def detuningscan():
        sim = RamanSimulation()
        average_energy = 3.0 * Er
        detunings = np.linspace(-1.5*sim.f_trap, 0.5*sim.f_trap, 20)
        print detunings.shape
        sample_idx = sim.numsteps - 1
        for delta in detunings:
            sim.detuning = delta
            lat = harmonicoscillator3d.harmonicoscillator3d(trap_frequencies, cutoffs)
            temperature = lat.temperature(average_energy)
            sim.init_populations[0:-1] = lat.populations_sum_over_all_but_first_frequency(temperature)
            sim.simulate_cooling(isPlot=False)
            
            axmean.plot(delta/1.0e3, sim.mean_phonon_number[sample_idx], marker = "o", linestyle = "None", color = 'b')
            print sim.mean_phonon_number[sample_idx]
            axmean.set_title("Mean phonon number vs. Detuning")
            axmean.set_ylabel("Mean phonon number")
            axmean.set_xlabel("Detuning from carrier (kHz)")
            
            axgnd.plot(delta/1.0e3, sim.populations[sample_idx, 0], marker = "o", linestyle = "None", color = 'b')
            axgnd.set_title("Ground state fraction vs. Detuning")
            axgnd.set_ylabel("Ground state fraction")
            axgnd.set_xlabel("Detuning from carrier (kHz)")
            
            axlost.plot(delta/1.0e3, sim.populations[sample_idx, sim._cutoff_state], marker = "o", linestyle = "None", color = 'b')
            axlost.set_title("Fraction lost vs. Detuning")
            axlost.set_ylabel("Fraction lost")
            axlost.set_xlabel("Detuning from carrier (kHz)")
            
    def durationscan():
        sim = RamanSimulation()
        average_energy = 3.0 * Er
        durations = np.linspace(0.05*sim.pi_time, 50*sim.pi_time, 20)
        sample_idx = sim.numsteps - 1
        for duration in durations:
            sim.duration = duration
            lat = harmonicoscillator3d.harmonicoscillator3d(trap_frequencies, cutoffs)
            temperature = lat.temperature(average_energy)
            sim.init_populations[0:-1] = lat.populations_sum_over_all_but_first_frequency(temperature)
            sim.simulate_cooling(isPlot=False)
            
            axmean.plot(duration*1.0e6, sim.mean_phonon_number[sample_idx], marker = "o", linestyle = "None", color = 'b')
            print sim.mean_phonon_number[sample_idx]
            axmean.set_title("duration vs. Detuning")
            axmean.set_ylabel("Mean phonon number")
            axmean.set_xlabel("Duration (us)")
            
            axgnd.plot(duration*1.0e6, sim.populations[sample_idx, 0], marker = "o", linestyle = "None", color = 'b')
            axgnd.set_title("Ground state fraction vs. Detuning")
            axgnd.set_ylabel("Ground state fraction")
            axgnd.set_xlabel("Duration (us)")
            
            axlost.plot(duration*1.0e6, sim.populations[sample_idx, sim._cutoff_state], marker = "o", linestyle = "None", color = 'b')
            axlost.set_title("Fraction lost vs. Duration")
            axlost.set_ylabel("Fraction lost")
            axlost.set_xlabel("Duration (us)")
    
    durationscan()
    plt.show()
    
    
    

