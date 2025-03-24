import numpy as np
import matplotlib.pyplot as plt

# Paramètres
L = 10.0        # Longueur de la région spatiale
T = 100.0       # Temps total de simulation
Nx = 100        # Nombre de points spatiaux
Nt = 500        # Nombre de pas de temps
k1 = 0.1        # Constante k1
D0 = 1e-3   # Coefficient de diffusion
C0 = 1.0        # Condition initiale de la concentration
dx = L / (Nx - 1)  # Pas d'espace
dt = T / Nt       # Pas de temps
alpha = D0 * dt / (dx ** 2)

# Fonction M1(x, t) (exemple sinus pour illustrer)
def M1(x, t):
    return 0.1 * (1 + 0.5 * np.sin(2 * np.pi * x / L) * np.cos(0.01 * t))

# Initialisation
X = np.linspace(0, L, Nx)
C = np.ones(Nx) * C0
C_new = np.zeros(Nx)
normes_L2 = []

# Matrice Tridiagonale (Implication diffusion)
A = np.eye(Nx) * (1 + 2 * alpha)
A += np.diag(-alpha * np.ones(Nx-1), k=1)
A += np.diag(-alpha * np.ones(Nx-1), k=-1)
A[0, :] = A[-1, :] = 0  # Conditions Dirichlet (bordures fixes)
A[0, 0] = A[-1, -1] = 1

# Boucle temporelle
for n in range(Nt):
    t = n * dt
    M = M1(X, t)
    B = C - k1 * M * C * dt  # Réaction explicite
    C_new = np.linalg.solve(A, B)  # Diffusion implicite
    C_new[0] = C_new[-1] = 0.0  # Condition Dirichlet
    C[:] = C_new[:]
    normes_L2.append(np.sqrt(np.sum(C**2)))

# Affichage
plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, T, Nt), normes_L2, label="Norme L2")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Temps (s)")
plt.ylabel("Norme L2")
plt.title("Évolution de la norme L2")
plt.grid(True)
plt.legend()
plt.show()
