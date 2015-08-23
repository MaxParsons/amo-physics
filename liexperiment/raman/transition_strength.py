'''
Created on Feb 20, 2015

@author: MP
'''
import amo.core.wigner
wigner = amo.core.wigner.Wigner
import numpy as np

def transition_through_single_excited(ii, jg0, jg1, je, fg0, fg1, fe, mg0, mg1, me, q0, q1):
    qs = [-1, 0, 1]
    answer = 0.0
    for idx in range(1,4):
        answer += q0[idx] * q1[idx] * (-1.0)**(fg0 + fg1 + 2.0 * fe + jg0 + je + 2.0 * ii - mg0 - me + 2) *\
        (2.0 * fe + 1) * np.sqrt((2 * fg0 + 1) * (2 * fg1 + 1)) *\
        wigner.wigner_3j(jg0, 1, je, -mg0, qs[idx], me) *\
        wigner.wigner_3j(je, 1, jg1, -me, qs[idx], mg1) *\
        wigner.wigner_6j(jg0, ii, fg0, fe, ii, je) *\
        wigner.wigner_6j(je, ii, fe, fg1, ii, jg1)
    return answer

wigner.wigner_6j(1, 0.5, 0.5, 0.5, 0.5, 0.5)

je = 1.5
fg0 = 0.5
fg1 = 1.5
fe = 1.5
mg0 = -0.5
mg1 = 0.5
me = 1.5
q0 = [0,0,1]
q1 = [0,1,0]
transition_through_single_excited(1.0, 0.5, 0.5, je, fg0, fg1, fe, mg0, mg1, me, q0, q1)