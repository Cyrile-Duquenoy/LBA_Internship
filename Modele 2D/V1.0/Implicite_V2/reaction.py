import numpy as np
import matplotlib.pyplot as plt
from laplacian_2D import *

def chemotaxis_flux(F, C, chi, dx, N):
    """
    Calcule ∇⋅(χ F ∇C) en différences finies avec condition de Neumann homogène.
    
    Paramètres :
    - F : Matrice (N, N) de concentration de fibroblastes
    - C : Matrice (N, N) de concentration de cytokines
    - chi : Coefficient de chimiotactisme (peut être un scalaire ou une matrice)
    - dx : Pas spatial
    - N : Taille du maillage (N x N)

    Retourne :
    - Matrice (N, N) représentant ∇⋅(χ F ∇C)
    """
    F = F.reshape(N, N)
    C = C.reshape(N, N)
    chi_F = chi * F  # Produit élément par élément si chi est une matrice

    # Allocation du résultat
    div_chiF_gradC = np.zeros((N, N))
    # Discrétisation par différences finies centrées
    for i in range(1, N-1):
        for j in range(1, N-1):
            dCdx = (C[i+1, j] - C[i-1, j]) / (2 * dx)
            dCdy = (C[i, j+1] - C[i, j-1]) / (2 * dx)

            d_chiF_dCdx = (chi_F[i+1, j] * dCdx - chi_F[i-1, j] * dCdx) / (2 * dx)
            d_chiF_dCdy = (chi_F[i, j+1] * dCdy - chi_F[i, j-1] * dCdy) / (2 * dx)

            div_chiF_gradC[i, j] = d_chiF_dCdx + d_chiF_dCdy

    # Conditions de Neumann homogènes (dérivée normale nulle sur les bords)
    div_chiF_gradC[0, :] = div_chiF_gradC[1, :]
    div_chiF_gradC[-1, :] = div_chiF_gradC[-2, :]
    div_chiF_gradC[:, 0] = div_chiF_gradC[:, 1]
    div_chiF_gradC[:, -1] = div_chiF_gradC[:, -2]
    return div_chiF_gradC.flatten()



def reaction(C,F, alpha, beta, gamma, delta, F_max, C_max, dx, N, chi):
    cytok_reac = alpha*C - beta*C
    
    #fib_reac = gamma*alpha*F*(1 - F/F_max) - delta*F - chemotaxis_flux(F, C, chi, dx, N)
    fib_reac = gamma*C*F/(F + F_max) - delta*F - chemotaxis_flux(F, C, chi, dx, N)
    
    return cytok_reac, fib_reac