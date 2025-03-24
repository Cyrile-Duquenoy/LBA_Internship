import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from module_perso import *

D=0.01 # Coeff. de diffusion
T=200 #Temps max en jours
L=1. # Largeur du carré
N=11
DT=50
dt=T/(DT-1)

x=np.linspace(0,L,N)
y=np.linspace(0,L,N)
X, Y = np.meshgrid(x, y)
dx=L/(N-1)
dy=dx

lbda=D*dt/(dx**2)

def f(x:float,y:float)->float:
    return x*(1-y)

Z = f(X, Y)

def Mat_x(N:int)->np.ndarray:
    size = N * N  # Taille de la matrice finale
    d = np.ones(size) * (1 - 4 * lbda)  # Diagonale principale
    d1 = np.ones(size - 1) * (-lbda)  # Diagonale adjacente
    
    # Supprimer les connexions entre les lignes différentes
    for i in range(1, N):
        d1[i * N - 1] = 0

    return np.diag(d) + np.diag(d1, -1) + np.diag(d1, 1)

def Mat_y(N:int)->np.ndarray:
    size = N * N  # Taille de la matrice finale
    d = np.ones(size - N) * (-lbda)  # Décalage entre les lignes
    return np.diag(d, -N) + np.diag(d, N)

def Mat_2D(N:int)->np.ndarray:
    return Mat_x(N) + Mat_y(N)

A=Mat_2D(3)
print(A)

# Résolution de l'équation en évoluant dans le temps
U = np.zeros((DT, N * N))  # Stocke les solutions à chaque instant
U[0, :] = Z  # Condition initiale

for t in range(1, DT):
    U[t, :] = A @ U[t - 1, :]  # Évolution temporelle

# Affichage de la solution finale
Z_final = U[-1, :].reshape(N, N)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z_final, cmap="viridis")

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Concentration")
ax.set_title("Évolution de la diffusion en 2D")

plt.show()

'''
# Création de la figure et de l'axe 3D
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

# Tracé de la surface
ax.plot_surface(X, Y, Z, cmap='viridis')

# Labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('f(X, Y)')
ax.set_title('Plot de f(x, y) = x(1 - y)')

plt.show()
'''