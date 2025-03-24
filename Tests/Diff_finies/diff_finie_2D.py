import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

L=1. # Largeur du carré
N=10 # Nombre de points de discrétisation par ligne
x=np.linspace(0, L, N)
dx=L/(N-1)
dy=dx
y=np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)
T=100 # Nombre total de pas de temps
dt=1 # Pas de temps
D=0.1 # Coefficient de diffusion

lbda=D*dt/(dx**2)

def u_init(x: float, y: float) -> float:
    if x == 0 or y == 0 or x == L or y == L:
        return 0
    else:
        return 1

def U0(x: np.ndarray) -> np.ndarray:
    U = np.zeros((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            U[i, j] = u_init(x[i], x[j])
    return U

def flux(x: float, y: float) -> float:
    return x * (1 - y)

def Flux(x: np.ndarray) -> np.ndarray:
    FF = np.zeros((len(x), len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            FF[i, j] = flux(x[i], x[j])
    return FF

def Mat_x(N: int) -> np.ndarray:
    d = np.ones(N*N) * (1 + 4 * lbda)
    d1 = np.ones((N*N)-1) * (-lbda)
    return np.diag(d) + np.diag(d1, -1) + np.diag(d1, 1)

def Mat_y(N: int) -> np.ndarray:
    d = np.ones(N*N-N) * lbda
    return np.diag(d, -N) + np.diag(d, N)

def Mat_2D(N: int) -> np.ndarray:
    return Mat_x(N) + Mat_y(N)

def solve_diffusion(N, T, dt):
    U = U0(x)
    M = Mat_2D(N)
    solutions = []
    
    for t in range(T+1):
        U = np.dot(M, U.flatten()).reshape(N, N) + Flux(x) * dt
        if t % 10 == 0:
            solutions.append(U.copy())
    
    return solutions

# Résolution
solutions = solve_diffusion(N, T, dt)

# Affichage
fig = plt.figure(figsize=(12, 8))
for idx, U in enumerate(solutions[:6]):
    ax = fig.add_subplot(2, 3, idx + 1, projection='3d')
    ax.plot_surface(X, Y, U, cmap="viridis")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Concentration")
    ax.set_title(f"t = {idx * 10}")

plt.tight_layout()
plt.show()
    
    


