'''
Created on Jul 27, 2015

@author: MP
'''
import numpy as np
import plastia.calculation.operator as op
import plastia.analysis.physics as phys
import matplotlib.pylab as plt


def ground_breit_rabi_plot():
    '''
    Construct ground electronic hyperfine Hamiltonian
    '''
    HhfGnd = op.Operator('F', [0, 0.5, 1])
    for J in np.arange(HhfGnd.momenta[0] + HhfGnd.momenta[1], -.1, -1):
        for F in np.arange(J + HhfGnd.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                if F == 1.5:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = phys.LithiumSix.mag_dipole_SOneHalf / 2  # HzGnd
                else:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = -phys.LithiumSix.mag_dipole_SOneHalf  # HzGnd
    
    '''
    Construct the ground electronic zeeman Hamiltonian (magnetic field factored out)
    '''
    
    HzGnd = op.Operator('L_S_I', [0, 0.5, 1])
    for mS in np.arange(-HzGnd.momenta[1], HzGnd.momenta[1] + .1, 1):
        for mI in np.arange(-HzGnd.momenta[2], HzGnd.momenta[2] + .1, 1):
            for mL in np.arange(-HzGnd.momenta[0], HzGnd.momenta[0] + .1, 1):
                HzGnd.representation[HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI), HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI)]\
     = phys.PhysicalConstants.gelectron * phys.PhysicalConstants.uB * mS + phys.PhysicalConstants.uB * mL + phys.PhysicalConstants.uN * mI
    HzGnd.ChangeOfBasis('F')
    
    B = np.arange(0, 528, 1)
    energiesGnd = np.zeros((HzGnd.length, len(B)))
    for n in np.arange(0, len(B), 1):
        w, v = np.linalg.eig(HhfGnd.representation + B[n] * HzGnd.representation)
        energiesGnd[:, n] = sorted(w)
    
    for k in np.arange(0, HzGnd.length, 1):
        plt.plot(B, energiesGnd[k, :])
    
    plt.show()
    
def excited_breit_rabi_plot():
    '''
    Construct hyperfine Hamiltonian
    '''
    Hhf = op.Operator('F', [1, 0.5, 1])
    for J in np.arange(Hhf.momenta[0] + Hhf.momenta[1], -.1, -1):
        for F in np.arange(J + Hhf.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                I = Hhf.momenta[2]
                C = F * (F + 1) - J * (J + 1) - I * (I + 1)
                if J == 0.5:
                    Hhf.representation[Hhf.FtoIndex(J, I, F, mF), Hhf.FtoIndex(J, I, F, mF)] = 0.5 * C * phys.LithiumSix.mag_dipole_POneHalf  # Hz
                else:
                    Hhf.representation[Hhf.FtoIndex(J, I, F, mF), Hhf.FtoIndex(J, I, F, mF)] = 0.5 * C * phys.LithiumSix.mag_dipole_PThreeHalf + \
                    0.375 * phys.LithiumSix.elec_quadrupole_PThreeHalf * C * (C + 1) / (I * (2 * I - 1) * J * (2 * J - 1))  # Hz
    
    '''
    Construct the fine structure Hamiltonian
    '''
    Hfs = op.Operator('J_I', [1, 0.5, 1])
    for J in np.arange(Hfs.momenta[0] + Hfs.momenta[1], -.1, -1):
        for mJ in np.arange(-J, J + .1, 1):
            for mI in np.arange(-Hfs.momenta[2], Hfs.momenta[2] + .1, 1):
                if J == 1.5:
                    Hfs.representation[Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI), Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI)] = phys.LithiumSix.fine_structure
                else:
                    Hfs.representation[Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI), Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI)] = 0.0
    Hfs.ChangeOfBasis('F')
    
    
    '''
    Construct the zeeman Hamiltonian (magnetic field factored out)
    '''
    
    Hz = op.Operator('L_S_I', [1, 0.5, 1])
    for mS in np.arange(-Hz.momenta[1], Hz.momenta[1] + .1, 1):
        for mI in np.arange(-Hz.momenta[2], Hz.momenta[2] + .1, 1):
            for mL in np.arange(-Hz.momenta[0], Hz.momenta[0] + .1, 1):
                Hz.representation[Hhf.L_S_ItoIndex(Hz.momenta[0], mL, Hz.momenta[1], mS, Hz.momenta[2], mI), Hhf.L_S_ItoIndex(Hz.momenta[0], mL, Hz.momenta[1], mS, Hz.momenta[2], mI)]\
     = phys.PhysicalConstants.gelectron * phys.PhysicalConstants.uB * mS + phys.PhysicalConstants.uB * mL + phys.PhysicalConstants.uN * mI
    Hz.ChangeOfBasis('F')
    
    B = np.arange(0, 2, .1)
    energies = np.zeros((Hz.length, len(B)))
    eigenvectors = np.zeros
    for n in np.arange(0, len(B), 1):
        w, v = np.linalg.eig(Hhf.representation + Hfs.representation + B[n] * Hz.representation)
        energies[:, n] = sorted(w)
    
    for k in np.arange(0, Hz.length, 1):
        plt.plot(B, energies[k, :])
    
    plt.show()
    
def imaging_transitions_high_field(field):
    for gstate in range(1, 4):
        freq = imaging_transition(field, gstate, 7) / 10**6
        print "|{}> --> |J=3/2, mJ=-3/2>: {} MHz".format(gstate, freq)
        
    for gstate in range(4, 7):
        freq = imaging_transition(field, gstate, 18) / 10**6
        print "|{}> --> |J=3/2, mJ=+3/2>: {} MHz".format(gstate, freq)
    

def print_microwave_transitions(field):
    '''
    Construct ground electronic hyperfine Hamiltonian
    '''
    HhfGnd = op.Operator('F', [0, 0.5, 1])
    for J in np.arange(HhfGnd.momenta[0] + HhfGnd.momenta[1], -.1, -1):
        for F in np.arange(J + HhfGnd.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                if F == 1.5:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = phys.LithiumSix.mag_dipole_SOneHalf / 2  # HzGnd
                else:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = -phys.LithiumSix.mag_dipole_SOneHalf  # HzGnd
    
    '''
    Construct the ground electronic zeeman Hamiltonian (magnetic field factored out)
    '''
    
    HzGnd = op.Operator('L_S_I', [0, 0.5, 1])
    for mS in np.arange(-HzGnd.momenta[1], HzGnd.momenta[1] + .1, 1):
        for mI in np.arange(-HzGnd.momenta[2], HzGnd.momenta[2] + .1, 1):
            for mL in np.arange(-HzGnd.momenta[0], HzGnd.momenta[0] + .1, 1):
                HzGnd.representation[HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI), HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI)]\
     = phys.PhysicalConstants.gelectron * phys.PhysicalConstants.uB * mS + phys.PhysicalConstants.uB * mL + phys.PhysicalConstants.uN * mI
    HzGnd.ChangeOfBasis('F')
    
    energiesGnd = np.zeros(HzGnd.length)
    
    w, v = np.linalg.eig(HhfGnd.representation + field * HzGnd.representation)
    energiesGnd[:] = sorted(w)
    
    for gstate in range(1,7):
        for estate in range (gstate + 1, 7):
            freq = (energiesGnd[estate - 1] - energiesGnd[gstate - 1])/10**6
            print "|{}> --> |{}>: {} MHz".format(gstate, estate, freq)

def imaging_transition(field, ground_state, excited_state):
    BEAT_FREQUENCY_OFFSET = 9887298340.0
    '''
    Construct ground electronic hyperfine Hamiltonian
    '''
    HhfGnd = op.Operator('F', [0, 0.5, 1])
    for J in np.arange(HhfGnd.momenta[0] + HhfGnd.momenta[1], -.1, -1):
        for F in np.arange(J + HhfGnd.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                if F == 1.5:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = phys.LithiumSix.mag_dipole_SOneHalf / 2  # HzGnd
                else:
                    HhfGnd.representation[HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF), HhfGnd.FtoIndex(J, HhfGnd.momenta[2], F, mF)] = -phys.LithiumSix.mag_dipole_SOneHalf  # HzGnd
    
    '''
    Construct the ground electronic zeeman Hamiltonian (magnetic field factored out)
    '''
    
    HzGnd = op.Operator('L_S_I', [0, 0.5, 1])
    for mS in np.arange(-HzGnd.momenta[1], HzGnd.momenta[1] + .1, 1):
        for mI in np.arange(-HzGnd.momenta[2], HzGnd.momenta[2] + .1, 1):
            for mL in np.arange(-HzGnd.momenta[0], HzGnd.momenta[0] + .1, 1):
                HzGnd.representation[HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI), HhfGnd.L_S_ItoIndex(HzGnd.momenta[0], mL, HzGnd.momenta[1], mS, HzGnd.momenta[2], mI)]\
     = phys.PhysicalConstants.gelectron * phys.PhysicalConstants.uB * mS + phys.PhysicalConstants.uB * mL + phys.PhysicalConstants.uN * mI
    HzGnd.ChangeOfBasis('F')
    
    energiesGnd = np.zeros(HzGnd.length)
    
    w, v = np.linalg.eig(HhfGnd.representation + field * HzGnd.representation)
    energiesGnd[:] = sorted(w)
    
    '''
    Construct Excited hyperfine Hamiltonian
    '''
    Hhf = op.Operator('F', [1, 0.5, 1])
    for J in np.arange(Hhf.momenta[0] + Hhf.momenta[1], -.1, -1):
        for F in np.arange(J + Hhf.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                I = Hhf.momenta[2]
                C = F * (F + 1) - J * (J + 1) - I * (I + 1)
                if J == 0.5:
                    Hhf.representation[Hhf.FtoIndex(J, I, F, mF), Hhf.FtoIndex(J, I, F, mF)] = 0.5 * C * phys.LithiumSix.mag_dipole_POneHalf  # Hz
                else:
                    Hhf.representation[Hhf.FtoIndex(J, I, F, mF), Hhf.FtoIndex(J, I, F, mF)] = 0.5 * C * phys.LithiumSix.mag_dipole_PThreeHalf + \
                    0.375 * phys.LithiumSix.elec_quadrupole_PThreeHalf * C * (C + 1) / (I * (2 * I - 1) * J * (2 * J - 1))  # Hz
    
    '''
    Construct the Excited fine structure Hamiltonian
    '''
    Hfs = op.Operator('J_I', [1, 0.5, 1])
    for J in np.arange(Hfs.momenta[0] + Hfs.momenta[1], -.1, -1):
        for mJ in np.arange(-J, J + .1, 1):
            for mI in np.arange(-Hfs.momenta[2], Hfs.momenta[2] + .1, 1):
                if J == 1.5:
                    Hfs.representation[Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI), Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI)] = phys.LithiumSix.fine_structure
                else:
                    Hfs.representation[Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI), Hfs.J_ItoIndex(J, mJ, Hfs.momenta[2], mI)] = 0.0
    Hfs.ChangeOfBasis('F')
    
    
    '''
    Construct the Excited zeeman Hamiltonian (magnetic field factored out)
    '''
    
    Hz = op.Operator('L_S_I', [1, 0.5, 1])
    for mS in np.arange(-Hz.momenta[1], Hz.momenta[1] + .1, 1):
        for mI in np.arange(-Hz.momenta[2], Hz.momenta[2] + .1, 1):
            for mL in np.arange(-Hz.momenta[0], Hz.momenta[0] + .1, 1):
                Hz.representation[Hhf.L_S_ItoIndex(Hz.momenta[0], mL, Hz.momenta[1], mS, Hz.momenta[2], mI), Hhf.L_S_ItoIndex(Hz.momenta[0], mL, Hz.momenta[1], mS, Hz.momenta[2], mI)]\
     = phys.PhysicalConstants.gelectron * phys.PhysicalConstants.uB * mS + phys.PhysicalConstants.uB * mL + phys.PhysicalConstants.uN * mI
    Hz.ChangeOfBasis('F')
    
    B = np.arange(0, 2, .1)
    energiesEx = np.zeros(Hz.length)
    w, v = np.linalg.eig(Hhf.representation + Hfs.representation + field * Hz.representation)
    energiesEx[:] = sorted(w)
    
    return energiesEx[excited_state - 1] - energiesGnd[ground_state - 1] - BEAT_FREQUENCY_OFFSET

def array_plot_test():
    a = np.array([1.0, 1.4, None, 3, 6, 2])
    pass
    
if __name__ == "__main__":
    #ground_breit_rabi_plot()
    #excited_breit_rabi_plot()
    energy = imaging_transitions_high_field(284)
    #print_microwave_transitions(10)
    