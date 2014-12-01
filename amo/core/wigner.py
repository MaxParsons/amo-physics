'''
Created on Nov 26, 2012

@author: Max
'''
"""wigner.py

Contains routines for computing Clebsch-Gordan coefficients, 3-j and
6-j symbols.

Note that these routines are *not* optimized in any way."""

from numpy import sqrt, array
from scipy.misc import factorial

class Wigner(object):

    # # CONSTANTS
    ##############
    
    # Limit on number of iterations for computing C-G coefficients
    CG_LIMIT = 50
    
    # Same but for the 6j calculation
    SIXJ_LIMIT = 50
    
    # # FUNCTIONS
    ##############
    
    @staticmethod
    def k_delta(x, y):
        """Kronecker delta function."""
        return (x == y)
    
    @staticmethod
    def triangle(a, b, c):
        """Compute the triangle coefficient Delta(a b c).
        See http://mathworld.wolfram.com/TriangleCoefficient.html."""
        return factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c) / factorial(a + b + c + 1)
    
    @classmethod
    def clebsch_gordan(cls, j1, j2, m1, m2, j, m):
        """Computes the Clebsch-Gordan coefficient
        <j1 j2; m1 m2|j1 j2; jm>.
    
        For reference see
        http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients."""
        kron = cls.k_delta(m, m1 + m2)
        if kron == 0:
            return kron
        c1 = sqrt((2 * j + 1) * factorial(j + j1 - j2) * factorial(j - j1 + j2) * \
                  factorial(j1 + j2 - j) / factorial(j1 + j2 + j + 1))
        c2 = sqrt(factorial(j + m) * factorial(j - m) * factorial(j1 - m1) * \
                  factorial(j1 + m1) * factorial(j2 - m2) * factorial(j2 + m2))
        c3 = 0.
        for k in range(cls.CG_LIMIT):
            use = True
            d = [0, 0, 0, 0, 0]
            d[0] = j1 + j2 - j - k
            d[1] = j1 - m1 - k
            d[2] = j2 + m2 - k
            d[3] = j - j2 + m1 + k
            d[4] = j - j1 - m2 + k
            prod = factorial(k)
            for arg in d:
                if arg < 0:
                    use = False
                    break
                prod *= factorial(arg)
            if use:
                c3 += (-1) ** k / prod
        return c1 * c2 * c3
    
    @classmethod
    def wigner_3j(cls, j1, j2, j3, m1, m2, m3):
        """Computes the Wigner 3-j symbol:
        ([ j1 j2 j3 ]
         [ m1 m2 m3 ])
        based on the appropriate Clebsch-Gordan coefficient."""
        return cls.clebsch_gordan(j1, j2, m1, m2, j3, -m3) * \
               (-1) ** (j1 - j2 - m3) / sqrt(2 * j3 + 1)
    
    @classmethod
    def wigner_6j(cls, j1, j2, j3, J1, J2, J3):
        """Computes the Wigner 6-j symbol:
        {[ j1 j2 j3]
         [ J1 J2 J3]}
        Notation follows Wolfram MathWorld."""
        def f(t):
            args = [t - j1 - j2 - j3, t - j1 - J2 - J3, t - J1 - j2 - J3,
                    t - J1 - J2 - j3, j1 + j2 + J1 + J2 - t, j2 + j3 + J2 + J3 - t,
                    j3 + j1 + J3 + J1 - t]
            args = factorial(array(args))
            result = args.prod()
            return result
                
        sixj = sqrt(cls.triangle(j1, j2, j3) * cls.triangle(j1, J2, J3) * \
                    cls.triangle(J1, j2, J3) * cls.triangle(J1, J2, j3))
        term = 0.
        for t in range(cls.SIXJ_LIMIT):
            g = f(t)
            if g <= 0:
                continue
            else:
                term += (-1) ** t * factorial(t + 1) / g
        sixj *= term
        return sixj
