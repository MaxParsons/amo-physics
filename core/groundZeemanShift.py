'''
Created on Mar 11, 2013

@author: Max

'''
import numpy as np
import plastia.calculation.operator as op
import plastia.analysis.physics as phys
import matplotlib.pylab as plt



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
energiesEx = np.zeros((HzGnd.length, len(B)))
for n in np.arange(0, len(B), 1):
    w, v = np.linalg.eig(HhfGnd.representation + B[n] * HzGnd.representation)
    energiesEx[:, n] = sorted(w)

for k in np.arange(0, HzGnd.length, 1):
    plt.plot(B, energiesEx[k, :])

plt.show()



