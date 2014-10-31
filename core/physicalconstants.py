'''
Created on Jun 3, 2014

@author: Max
'''
import numpy as np

class PhysicalConstants(object):
    h = 6.626068e-34# m**2 * kg / s
    hbar = 1.05457148e-34# m**2 * kg / s
    amu = 1.66053886e-27# kg
    e = 1.60217646e-19# C
    kb = 1.3806503e-23# m**2 * kg / s**2 / K
    c = 299792458# m**2 / s
    g = 9.80559# m / s**2 in Cambridge
    uB = 1.3996245e10# Hz/T (Bohr Magneton)
    gelectron = 2.002319304# electron g-factor
    uN = 7.62259277e6# Hz/T (Nuclear Magneton)
    epsilon0 = 8.854187817e-12# farads/meter vacuum permittivity


class LithiumSix(object):
    mass = 9.96323318e-27# kg

    wavelength_d1 = 670.992421e-9# m
    wavelength_d2 = 670.977338e-9# m

    frequency_d1 = 446.789634e12# Hz
    frequency_d2 = 446.799677e12# Hz
    fine_structure = 10.053044e9# Hz

    linewidth = 5.8724e6# Hz, not s**(-1)!

    mag_dipole_SOneHalf = 152.1368407e6# Hz, magnetic dipole hyperfine constant for ground state
    mag_dipole_POneHalf = 17.386e6# Hz, magnetic dipole hyperfine constant for 2P_1/2 state
    mag_dipole_PThreeHalf = -1.155e6# Hz, magnetic dipole hyperfine constant for 2P_3/2
    elec_quadrupole_PThreeHalf = -0.10e6# Hz, magnetic dipole hyperfine constant for 2P_3/2 state
    
    def polarizability(self, wavelength, state):# FIXME: Polarizabilities for lithium
        pass
    
    def scattering_length(self, field, state0, state1):# FIXME: Scattering lengths for lithium
        pass

class Rubidium87(object):
    mass = 1.44316077e-25# kg
    
