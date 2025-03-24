import numpy as np

##################################
######### CONSTANTES #############
##################################

# Paramètres
L = 10.0        # Longueur de la région spatiale
T = 100.0       # Temps total de simulation
Nx = 100        # Nombre de points spatiaux
Nt = 500        # Nombre de pas de temps
k1 = 0.1        # Constante k1
M1 = 1.0        # Constante M1
D0 = 0.01       # Coefficient de diffusion
C0 = 1.0        # Condition initiale de la concentration