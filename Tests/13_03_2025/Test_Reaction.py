import numpy as np
import matplotlib.pyplot

L=1
N=50

x=np.linspace(0,L,N)


#Source
def E(U):
    return True

#Production
def P(U):
    return True

#Apoptose
def D(U):
    return True

# Terme de Réaction
def R(U):
    return E(U)+P(U)+D(U)

