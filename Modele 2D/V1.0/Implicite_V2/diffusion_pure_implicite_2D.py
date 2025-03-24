import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal

from reaction_term import *
from param import *
from neumann import *

# Discret Param
N = 30
M = 9
L = 1.0
T = 200
dt = 0.001
t = np.arange(0, T + dt, dt)


x = np.linspace(0, L, N)
y = np.linspace(0, L, M)
X, Y = np.meshgrid(x, y)

dx = L / (N - 1)

# Param
rho = 1.2 * 10e-2
rad = N // 2
x0, y0 = N // 2, N // 2

D = 1.47*10e-4
lbda = D * dt / dx**2

# Initial Fibrosis
def u_init(N):
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad: 
                U[i, j] = rho
    return U

# Sparse_Matrix
def Mat_sparse(N):
    size = N * N
    diag_main = np.ones(size) * (1 + 4 * lbda)
    diag_x = np.ones(size - 1) * -lbda
    diag_y = np.ones(size - N) * -lbda
    

    for i in range(1, N):
        diag_x[i * N - 1] = 0  # Gestion des bords

    A = sp.diags([diag_main, diag_x, diag_x, diag_y, diag_y], [0, -1, 1, -N, N], format='csr')

    # Neumann Conditions
    for i in range(N):  
        for j in range(N):
            index = i * N + j  # Index global dans le vecteur aplati
            
            if i == 0:  # Bord bas
                A[index, index + N] += lbda  # Copie du point adjacent
            if i == N - 1:  # Bord haut
                A[index, index - N] += lbda
            if j == 0:  # Bord gauche
                A[index, index + 1] += lbda
            if j == N - 1:  # Bord droit
                A[index, index - 1] = lbda
    
    return A

# Initialisation
U = u_init(N).flatten()
A = Mat_sparse(N)

l2_norms = []



# Évolution temporelle
for n in range(len(t)):
    U = spsolve(A,U)  # Produit matrice creuse/vecteur
    l2_norms.append(np.linalg.norm(U, ord=1))

    if n % 100 == 0:
        plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
        plt.colorbar(label="Concentration de fibroblastes")
        plt.title(f"Évolution à t = {t[n]:.2f} j")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.pause(0.1)
plt.show()

# Plot de la norme L2
plt.figure()
plt.plot(t[:len(l2_norms)], l2_norms, label="Quantité totale")
plt.xlabel("Temps (j)")
plt.ylabel("Norme L2")
#plt.yscale('log')
#plt.xscale('log')
plt.title("Diffusion pure")
plt.legend()
plt.grid()
plt.show()