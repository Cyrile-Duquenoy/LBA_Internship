import numpy as np
from param import *

'''
def P(u, aplha, u_max):
    return alpha*u*(1 - u/u_max)

def D(u,beta):
    return beta*u

def R(u):
    return P(u) - D(u)
'''

def P(U):
    I1 = (T_b/(K_Tb + T_b) + I13/(K_I13 + I13))
    return lbda_fE*I1*(E / (K_E + E))*U

def D(U):
    return d_f*U

def fib_to_myo(U):
    I1 = T_b/(K_Tb + T_b)
    I2 = G/(K_G+G)
    return (lbda_mfT*I1 + lbda_mfG*I2)*U

def P2(U):
    I = H_A/(K_HA + H_A)
    return lbda_fHA*I*U

def source(U):
    return lbda_Ef*E0

def R(U):
    return source(U) + P(U) + P2(U) - fib_to_myo(U) - D(U)