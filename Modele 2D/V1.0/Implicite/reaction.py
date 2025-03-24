import numpy as np
import matplotlib.pyplot as plt


def reaction(C,F, alpha, beta, gamma, delta, F_max, C_max, dx):
    Xsi=2
    cytok_reac = 0
    fib_reac = gamma*C*(F/(F + F_max))  + F*np.gradient(Xsi*np.gradient(C)) - delta*F
    return cytok_reac, fib_reac