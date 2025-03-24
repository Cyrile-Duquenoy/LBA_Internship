import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal

from reaction_term import *
from param import *

from laplacian_2D import *
from init import *
from reaction import *


# Discret Param
N = 30
M=N
L = 1.0
T = 50
dt = 0.01
t = np.arange(0, T + dt, dt)

x = np.linspace(0, L, N)
y = np.linspace(0, L, M)
X, Y = np.meshgrid(x, y)

dx = L / (N - 1)

# Param
rho = 1.2e-2
rad = N // 2
x0, y0 = N // 2, N // 2

D_F = 1.47e-3 # Fibro. Diffusion. Coeff.
D_C = 1.08e-2 # Cytok. Diffusion Coeff.


lbda_F = D_F * dt / dx**2
lbda_C = D_C * dt / dx**2

alpha = 1.39e-2 #Cyto. Prod.
beta = 55.45e-3 #Cyto. Death.
gamma = lbda_fE #Fibro. Prod.
delta = d_f #Apoptose Coeff.

A_F = Laplacian_2d(N, lbda_F)
A_C = Laplacian_2d(N, lbda_C)



size=N**2

if __name__ == "__main__":

    C0, F0 = initialisation(N, rho, L, rad)
    
    cytok_norm = []
    fib_norm = []
    # Évolution temporelle
    for n in range(1,len(t)):

        F0 = spsolve(A_F,F0)  # Produit matrice creuse/vecteur
        C0 = spsolve(A_C,C0)

        # Reaction Term
        R_C, R_f = reaction(C0, F0, alpha, beta, gamma, delta, rho/2, rho/2, dx)
        
        
        C0+=dt*R_C
        F0+=dt*R_f
        
        cytok_norm.append(np.linalg.norm(C0, ord=1)/size)
        fib_norm.append(np.linalg.norm(F0, ord=1)/size)
        
        

        if n % 100 == 0:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 ligne, 2 colonnes
            im1 = axes[0].imshow(C0.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
            axes[0].set_title(f"Cytokines à t = {t[n]:.2f} j")
            axes[0].set_xlabel("x")
            axes[0].set_ylabel("y")
            fig.colorbar(im1, ax=axes[0], label="Concentration de cytokines")

            im2 = axes[1].imshow(F0.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
            axes[1].set_title(f"FIbroblastes à t = {t[n]:.2f} j")
            axes[1].set_xlabel("x")
            axes[1].set_ylabel("y")
            fig.colorbar(im2, ax=axes[1], label="Concentration de fibroblastes")

            plt.tight_layout()  # Ajuste l'espacement pour éviter le chevauchement
            plt.show()
    plt.show()
    
    ##############################################################
    plt.figure()
    plt.plot(t[:len(cytok_norm)], cytok_norm, label="\u03C3 (Cytokine)/Jour")
    plt.xlabel("Temps (j)")
    plt.ylabel("Cytokine moyenne / Unité")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.title("Cytokine")
    plt.legend()
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.plot(t[:len(fib_norm)], fib_norm, label="\u03C3 (Fibroblaste)/Jour")
    plt.xlabel("Temps (j)")
    plt.ylabel("Fibroblaste moyenne / Unité")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.title("Fibroblaste")
    plt.legend()
    plt.grid()
    plt.show()