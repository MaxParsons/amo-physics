'''
Created on Feb 19, 2015

@author: Max
'''
import liexperiment.raman.coherent_population_transfer as cpt
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
from itertools import product

def export_figure_numerical_index(filename, fig):
    head, tail = os.path.split(filename)
    fig_nums = [int(fname[-8:-4]) for fname in os.listdir(head) if fname.split('_', 1)[0] == tail]
    if not fig_nums:
        next_num = 0
    else:
        next_num = np.max(np.array(fig_nums)) + 1
    newname = tail + "_" + "{:0>4d}".format(next_num)
    fig.savefig(os.path.join(head, newname + ".svg"))
    
def spectrum_constant_pulse():
    fig_directory = "C:\\Users\\Max\\amo-physics\\liexperiment\\raman\\coherent_population_transfer\\constant_detuning_rabi"
    subname = "spectrum"
    
    raman = cpt.RamanTransition()
    detunings = np.linspace(-2.0e6, 2.0e6, 100)
    four_pops = np.zeros_like(detunings)
    nbars = np.zeros_like(detunings)
    
    raman.n_vibrational = 5;
    raman.initial_state = np.zeros(2 * raman.n_vibrational, dtype="complex64")
    raman.constant_rabi = 300.0e3
    raman.anharmonicity = 26.0e3
    raman.simulation_duration = 10.0e-6
    raman.simulation_nsteps = 50
    raman.trap_frequency = 1.0e6
    raman.lamb_dicke = 0.28
    raman.initial_state[0] = np.sqrt(0.7)
    raman.initial_state[1] = np.sqrt(0.3)
    
    fig, ax = plt.subplots(1, 1)
    fig.name = "spectrum"
    ax.set_title("simulated raman spectrum\n ")
    ax.set_xlabel("detuning (kHz)")
    ax.set_ylabel("population in |4> (blue)\n nbar (black)")
    
    for idx, detuning in enumerate(detunings):
        print "idx = " + str(idx)
        raman.constant_detuning = detuning
        raman.compute_dynamics()
        four_pops[idx] = raman.pops_excited[-1]
        nbars[idx] = raman.nbars[-1]
    ax.plot(detunings / 1.0e3, four_pops, color="b", marker="o")
    ax.plot(detunings / 1.0e3, nbars, color="k", marker="o")
    export_figure_numerical_index(os.path.join(fig_directory, fig.name), fig)
    
    plt.show()
    
def rabi_flopping():
    fig_directory = "C:\\Users\\Max\\amo-physics\\liexperiment\\raman\\coherent_population_transfer\\constant_detuning_rabi"
    subname = "spectrum"
    
    raman = cpt.RamanTransition()
    raman.constant_detuning = -1.00e6
    
    raman.n_vibrational = 5;
    raman.initial_state = np.zeros(2 * raman.n_vibrational, dtype="complex64")
    raman.constant_rabi = 100.0e3
    raman.anharmonicity = 0.0e3
    raman.simulation_duration = 100.0e-6
    raman.simulation_nsteps = 100
    raman.trap_frequency = 1.0e6
    raman.lamb_dicke = 0.28
    raman.initial_state[0] = np.sqrt(0.7)
    raman.initial_state[1] = np.sqrt(0.3)
    raman.compute_dynamics()
    
    fig, ax = plt.subplots(1, 1)
    ax.set_title("populations")
    ax.set_xlabel("time")
    ax.set_ylabel("populations")
    plt.plot(raman.times, raman.pops_excited)
    plt.show()
    
def test():
    fig_directory = "C:\\Users\\Max\\amo-physics\\liexperiment\\raman\\coherent_population_transfer\\constant_detuning_rabi"
    subname = "spectrum"
    
    raman = cpt.RamanTransition()
    detunings = np.linspace(-2.0e6, 2.0e6, 30)
    four_pops = np.zeros_like(detunings)
    nbars = np.zeros_like(detunings)
    
    raman.n_vibrational = 3;
    raman.initial_state = np.zeros(2 * raman.n_vibrational, dtype="complex64")
    raman.constant_rabi = 100.0e3
    raman.anharmonicity = 26.0e3
    raman.simulation_duration = 10.0e-6
    raman.simulation_nsteps = 50
    raman.trap_frequency = 1.0e6
    raman.lamb_dicke = 0.28
    raman.initial_state[0] = np.sqrt(1.0)
    raman.initial_state[1] = np.sqrt(0.0)
    
    fig, ax = plt.subplots(1, 1)
    fig.name = "spectrum"
    ax.set_title("simulated raman spectrum\n ")
    ax.set_xlabel("detuning (kHz)")
    ax.set_ylabel("population in |4> (blue)\n nbar (black)")
    
    raman.constant_detuning = 1.0e6
    raman.compute_quantum_numbers()
    print raman.hamiltonian(2.2e-6)

if __name__ == "__main__":
    # test()
    spectrum_constant_pulse()
    # rabi_flopping()
