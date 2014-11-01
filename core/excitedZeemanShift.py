'''
Created on Mar 11, 2013

@author: Max

'''
import numpy as np
from atoms.hilberspace import Operator
import core.physicalconstants as phys
import matplotlib.pylab as plt

'''
Construct hyperfine Hamiltonian
'''
Hhf = Operator('F', [1, 0.5, 1])
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
Hfs = Operator('J_I', [1, 0.5, 1])
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

Hz = Operator('L_S_I', [1, 0.5, 1])
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



