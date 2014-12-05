'''
Created on Dec 5, 2014

@author: Max
'''
import numpy as np
import matplotlib.pyplot as plt
from amo.core import odesolver

if __name__ == "__main__":
    # SI units everywhere
    # experiment parameters
    f_trap = 0.65e6# trap frequency
    cutoff_state = 8.0# first state to tunnel, indexed from 0
    pi_time = 2.0e-6# raman carrier pi-time
    pump_time = 10.0e-6# pumping time
    photons_per_raman = 3.0# recoil heating for fudge factor
    
    # constants
    F_RECOIL = 75.0e3
    F_ANHARMONICITY = 26.0e3 
    
    # calculated from experiment parameters
    eta_0 = np.sqrt(F_RECOIL / f_trap)
    eta = 0.19 * np.sqrt(1.0e6 / f_trap)
    eps = photons_per_raman * eta_0 ** 2
    
    # SIMULATION
    # setup arrays
    init_populations = np.zeros((cutoff_state + 1,))
    init_populations[0] = 0.5
    init_populations[2] = 0.5
    transition_matrix = np.zeros((cutoff_state + 1, cutoff_state + 1))
    
    def lineshape(detuning, raman_rabi, decay_rate):
        s0 = 2.0 * raman_rabi ** 2 / decay_rate ** 2
        return decay_rate / 2.0 * s0 / (1 + s0 + (2.0 * detuning / decay_rate) ** 2)
    
    def sideband_frequencies(n0, nf):
        if n0 < nf:
            return f_trap * (nf - n0) - np.sign(nf - n0) * F_ANHARMONICITY * np.sum(np.arange(n0, nf))
        else:
            return f_trap * (nf - n0) - np.sign(nf - n0) * F_ANHARMONICITY * np.sum(np.arange(nf, n0))
    
    def transition_rate_neighbor_only(n0, nrate, detuning):
        decay_rate = 1.0 / pump_time
        if nrate == cutoff_state:
            if (nrate - n0) == 1:
                up_rabi = np.sqrt(n0 + 1) * eta * np.pi / pi_time
                stay_rabi = np.pi / pi_time
                f_heat = sideband_frequencies(n0, n0 + 1)
                return (1 + eps) * (lineshape(detuning - f_heat, up_rabi, decay_rate) + eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
            elif nrate == n0:
                down_rabi = np.sqrt(n0) * eta * np.pi / pi_time
                up_rabi = np.sqrt(n0 + 1) * eta * np.pi / pi_time
                stay_rabi = np.pi / pi_time
                f_heat = sideband_frequencies(n0, n0 + 1)
                f_cool = sideband_frequencies(n0, n0 - 1)
                return (1 - eps) * (-lineshape(detuning - f_cool, down_rabi, decay_rate) - eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
            else:
                return 0.0
        else:
            if (nrate - n0) == 1:
                up_rabi = np.sqrt(n0 + 1) * eta * np.pi / pi_time
                stay_rabi = np.pi / pi_time
                f_heat = sideband_frequencies(n0 , n0 + 1)
                return (1 + eps) * (lineshape(detuning - f_heat, up_rabi, decay_rate) + eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
            elif (nrate - n0) == -1:
                down_rabi = np.sqrt(n0) * eta * np.pi / pi_time
                stay_rabi = np.pi / pi_time
                f_cool = sideband_frequencies(n0, n0 - 1)
                return (1 - eps) * (lineshape(detuning - f_cool, down_rabi, decay_rate) + eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
            elif nrate == n0:
                down_rabi = np.sqrt(n0) * eta * np.pi / pi_time
                up_rabi = np.sqrt(n0 + 1) * eta * np.pi / pi_time
                stay_rabi = np.pi / pi_time
                f_heat = sideband_frequencies(n0, n0 + 1)
                f_cool = sideband_frequencies(n0, n0 - 1)
                if n0 == 0:
                    return (1 + eps) * (-lineshape(detuning - f_heat, up_rabi, decay_rate) - eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
                else:
                    return (1 + eps) * (-lineshape(detuning - f_heat, up_rabi, decay_rate) - eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate)) + \
                        (1 - eps) * (-lineshape(detuning - f_cool, down_rabi, decay_rate) - eta_0 ** 2 * lineshape(detuning, stay_rabi, decay_rate))
            else:
                return 0.0
    
    for n0 in np.arange(0, len(init_populations)):
        for nrate in np.arange(0, len(init_populations)):
            transition_matrix[nrate, n0] = transition_rate_neighbor_only(n0, nrate, -f_trap)

    print np.sum(transition_matrix, 0)
    r = odesolver.rateequation(transition_matrix, init_populations)
    r.solve(1.0e-3, numsteps=50)
    r.result.plot()
    plt.show()
    
    

