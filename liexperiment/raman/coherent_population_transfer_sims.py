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
        print np.array(fig_nums)
        next_num = np.max(np.array(fig_nums)) + 1
    newname = tail + "_" + "{:0>4d}".format(next_num)
    print newname
    fig.savefig(os.path.join(head, newname + ".svg"))
    
def spectrum_constant_pulse():
    fig_directory = "C:\\Users\\Max\\amo-physics\\liexperiment\\raman\\coherent_population_transfer\\constant_detuning_rabi"
    subname = "spectrum"
    
    raman = cpt.RamanTransition()
    detunings = np.linspace(-2.0e6, 2.0e6, 100)
    four_pops = np.zeros_like(detunings)
    nbars = np.zeros_like(detunings)
    
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

if __name__ == "__main__":
    spectrum_constant_pulse()
    
