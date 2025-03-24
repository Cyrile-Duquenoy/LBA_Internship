import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Pour le tracé 3D

# Paramètres
N = 50  # Taille de la grille (ajuste si nécessaire)
L = 1    # Longueur de l'espace
T = 0.01  # Temps total réduit pour test rapide
dx = L / (N - 1)
dt = 0.001

D = 0.001  # Coefficient de diffusion
lbda = D * dt / dx**2

rho = 1.2 * 10e-2  # Densité initiale
rad = N // 2
x0, y0 = N // 2, N // 2

# Initialisation de U (matrice aplatie)
def u_init(N):
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad: 
                U[i, j] = rho
    return U

U = u_init(N).flatten()

# Construction de la matrice creuse pour la diffusion
def Mat_sparse(N):
    size = N * N
    diag_main = np.ones(size) * (1 - 4 * lbda)
    diag_x = np.ones(size - 1) * lbda
    diag_y = np.ones(size - N) * lbda
    
    for i in range(1, N):
        diag_x[i * N - 1] = 0  # Gestion des bords
    
    A = sp.diags([diag_main, diag_x, diag_x, diag_y, diag_y], [0, -1, 1, -N, N], format='csr')
    
    return A

A = Mat_sparse(N)

# Grille pour le tracé 3D
X = np.linspace(0, L, N)
Y = np.linspace(0, L, N)
X, Y = np.meshgrid(X, Y)

# Initialisation du tracé
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Évolution temporelle avec affichage dynamique
t = np.arange(0, T + dt, dt)
for n in range(len(t)):
    U = A @ U  # Produit matrice creuse/vecteur

    if n % 100 == 0:  # Mise à jour toutes les 1000 itérations pour éviter un affichage trop lent
        #ax.cla()  # Efface l'ancienne figure
        U_matrix = U.reshape(N, N)

        surf = ax.plot_surface(X, Y, U_matrix, cmap='viridis', edgecolor='k')
        ax.set_title(f"Évolution à t = {t[n]:.5f} s")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("Concentration")

        #plt.pause(0.1)  # Pause pour l'animation

plt.show()



