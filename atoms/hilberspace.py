'''
Created on Mar 11, 2013

@author: Max

TODO: Change basis definitions to enum or add exceptions for giberish labels
'''

import numpy as np
from core.wigner import Wigner

class Operator(object):
    """
    Represents an operator in the Hilbert space of a single atomic electronic manifold.
    """

    def __init__(self, basis, momenta):
        """
        create a new state representation

        @param representation: A matrix representation of some operator.  The ordering of the array is from most negative to
        most positive angular momentum with the order of different angular momentum components the same as the 'basis'
        labels.  e.g. for J=3/2, I=1/2 (J,mJ,I,mI): vector =
        [(1/2,-1/2,1/2,-1/2),(1/2,-1/2,1/2,1/2),(1/2,1/2,1/2,-1/2),(1/2,1/2,1/2,1/2),...,
        (3/2,-3/2,1/2,-1/2),(3/2,-3/2,1/2,1/2),...

        @param basis: 'F','J_I','L_S_I'

        @param momenta: Maximum possible L,S,I.
        """
        assert (basis in ('F', 'J_I', 'L_S_I')), "Not a valid basis!"

        self.length = (2 * momenta[0] + 1) * (2 * momenta[1] + 1) * (2 * momenta[2] + 1)
        self.representation = np.zeros((self.length, self.length));
        self.basis = basis;
        self.momenta = momenta;

    def ChangeOfBasis(self, newBasis):
        """
        change to a different basis

        @param newBasis: 'F','J_I','L_S_I'
        """

        JToF = np.zeros((self.length, self.length));
        for J in np.arange(self.momenta[0] + self.momenta[1], -.1, -1):
            for mJ in np.arange(-J, J + .1, 1):
                for F in np.arange(J + self.momenta[2], -.1, -1):
                    for mF in np.arange(-F, F + .1, 1):
                        for mI in np.arange(-self.momenta[2], self.momenta[2] + .1, 1):
                            JToF[Operator.FtoIndex(J, self.momenta[2], F, mF), Operator.J_ItoIndex(J, mJ, self.momenta[2], mI)]\
 = Wigner.clebsch_gordan(J, self.momenta[2], mJ, mI, F, mF);

        LToJ = np.zeros((self.length, self.length));
        for J in np.arange(self.momenta[0] + self.momenta[1], -.1, -1):
            for mJ in np.arange(-J, J + .1, 1):
                for mS in np.arange(-self.momenta[1], self.momenta[1] + .1, 1):
                    for mI in np.arange(-self.momenta[2], self.momenta[2] + .1, 1):
                        for mL in np.arange(-self.momenta[0], self.momenta[0] + .1, 1):
                            LToJ[Operator.J_ItoIndex(J, mJ, self.momenta[2], mI), Operator.L_S_ItoIndex(self.momenta[0],
 mL, self.momenta[1], mS, self.momenta[2], mI)] = Wigner.clebsch_gordan(self.momenta[0], self.momenta[1], mL, mS, J, mJ);

        if self.basis == 'L_S_I' and newBasis == 'J_I':
            self.representation = np.dot(LToJ, np.dot(self.representation, np.linalg.inv(LToJ)))

        if self.basis == 'J_I' and newBasis == 'F':
            self.representation = np.dot(JToF, np.dot(self.representation, np.linalg.inv(JToF)))

        if self.basis == 'L_S_I' and newBasis == 'F':
            self.representation = np.dot(np.dot(JToF, LToJ), np.dot(self.representation, np.linalg.inv(np.dot(JToF, LToJ))))

        if self.basis == 'J_I' and newBasis == 'L_S_I':
            self.representation = np.dot(np.linalg.inv(LToJ), np.dot(self.representation, LToJ))

        if self.basis == 'F' and newBasis == 'J_I':
            self.representation = np.dot(np.linalg.inv(JToF), np.dot(self.representation, JToF))

        if self.basis == 'F' and newBasis == 'L_S_I':
            self.representation = np.dot(np.linalg.inv(np.dot(JToF, LToJ)), np.dot(self.representation, np.dot(JToF, LToJ)))

        self.basis = newBasis

    @staticmethod
    def J_ItoIndex(J, mJ, I, mI):
        return (2 * I + 1) * np.sum(2 * (np.arange(J - 1, 0, -1)) + 1) + (2 * I + 1) * (J + mJ) + I + mI;

    @staticmethod
    def FtoIndex(J, I, F, mF):
        j = J
        lowerManifoldStates = 0;
        while j - 1 > 0:
            lowerManifoldStates += np.sum(2 * (np.arange(j + I - 1, 0, -1)) + 1)
            j = j - 1
        
        return np.sum(2 * (np.arange(F - 1, 0, -1)) + 1) + F + mF + lowerManifoldStates

    @staticmethod
    def L_S_ItoIndex(L, mL, S, mS, I, mI):
        return (2 * I + 1) * (2 * S + 1) * (L + mL) + (2 * I + 1) * (S + mS) + I + mI;

class StateVector(object):
    """
    Represents an atomic or molecular angular momentum state vector.
    """
    
    def __init__(self, representation, basis, momenta):
        """
        create a new state representation
        
        @param representation: An array that gives the state vector.  The ordering of the array is from most negative to 
        most positive angular momentum with the order of different angular momentum components the same as the 'basis'
        labels.  e.g. for J=3/2, I=1/2 (J,mJ,I,mI): vector = 
        [(1/2,-1/2,1/2,-1/2),(1/2,-1/2,1/2,1/2),(1/2,1/2,1/2,-1/2),(1/2,1/2,1/2,1/2),...,
        (3/2,-3/2,1/2,-1/2),(3/2,-3/2,1/2,1/2),...
        
        @param basis: 'F','J_I','L_S_I'
        
        @param momenta: Maximum possible L,S,I.
        """
        #EXCEPTIONS!!!
        
        self.representation = representation;
        self.basis = basis;
        self.momenta = momenta;
     
    def ChangeOfBasis(self, newBasis):
        """
        change to a different basis
        
        @param newBasis: 'F','J_I','L_S_I'
        """

        JToF = np.zeros(2 * self.length + 1, 2 * self.length + 1);
        for F in np.arange(self.momenta[0] + self.momenta[1] + self.momenta[2], -.1, -1):
            for mF in np.arange(-F, F + .1, 1):
                for J in np.arange(self.momenta[0] + self.momenta[1], -.1, -1):
                    for mJ in np.arange(-J, J + .1, 1):
                        for mI in np.arange(-self.momenta[2], self.momenta[2] + .1, 1):
                            JToF[StateVector.FtoIndex(J, self.momenta[2], F, mF), StateVector.J_ItoIndex(J, mJ, self.momenta[2], mI)]\
 = Wigner.clebsch_gordan(J, self.momenta[2], mJ, mI, F, mF);
 
        LToJ = np.zeros(2 * self.length + 1, 2 * self.length + 1);
        for J in np.arange(self.momenta[0] + self.momenta[1] + self.momenta[2], -.1, -1):
            for mJ in np.arange(-J, J + .1, 1):
                for mS in np.arange(-self.momenta[1], self.momenta[1] + .1, 1):
                    for mI in np.arange(-self.momenta[2], self.momenta[2] + .1, 1):
                        for mL in np.arange(-self.momenta[0], self.momenta[0] + .1, 1):
                                JToF[StateVector.J_ItoIndex(J, mJ, self.momenta[1], mI), StateVector.L_S_ItoIndex(self.momenta[0],
mL, self.momenta[1], mS, self.momenta[2], mI)] = Wigner.clebsch_gordan(self.momenta[0], self.momenta[1], mL, mS, J, mJ);

        if self.basis == 'L_S_I' and newBasis == 'J_I':
            self.representation = np.dot(LToJ, self.representation)
        
        if self.basis == 'J_I' and newBasis == 'F':
            self.representation = np.dot(JToF, self.representation)
        
        if self.basis == 'L_S_I' and newBasis == 'F':
            self.representation = np.dot(LToJ, np.dot(JToF, self.representation))
        
        if self.basis == 'J_I' and newBasis == 'L_S_I':
            self.representation = np.dot(np.linalg.inv(LToJ), self.representation)
        
        if self.basis == 'F' and newBasis == 'J_I':
            self.representation = np.dot(np.linalg.inv(JToF), self.representation)
        
        if self.basis == 'F' and newBasis == 'L_S_I':
            self.representation = np.dot(np.linalg.inv(LToJ), np.dot(np.linalg.inv(JToF), self.representation))
        
        self.basis = newBasis;
                        
    @staticmethod
    def J_ItoIndex(J, mJ, I, mI):
        return (2 * I + 1) * np.sum(2 * (np.arange(J - 1, 0, -1)) + 1) + (2 * I + 1) * (J + mJ) + I + mI;
    
    @staticmethod
    def FtoIndex(J, I, F, mF):
        j = J
        lowerManifoldStates = 0;
        while j - 1 > 0:
            lowerManifoldStates += np.sum(2 * (np.arange(j + I - 1, 0, -1)) + 1)
            j = j - 1
        
        return np.sum(2 * (np.arange(F - 1, 0, -1)) + 1) + F + mF + lowerManifoldStates
        
    @staticmethod
    def L_S_ItoIndex(L, mL, S, mS, I, mI):
        return (2 * I + 1) * (2 * S + 1) * (L + mL) + (2 * I + 1) * (S + mS) + I + mI;
                        
