import numpy as np
import matplotlib.pyplot as plt

# Paramètres
Lx, Ly = 10.0, 10.0   # Longueur de la région spatiale en x et y
T = 10.0             # Temps total de simulation
Nx, Ny = 50, 50       # Nombre de points spatiaux en x et y
Nt = 500              # Nombre de pas de temps
k1 = 0.1              # Constante k1
D0 = 8.64e-7          # Coefficient de diffusion
C0 = 1.0              # Condition initiale de la concentration
dx, dy = Lx / (Nx - 1), Ly / (Ny - 1)  # Pas d'espace en x et y
dt = T / Nt           # Pas de temps
alpha_x = D0 * dt / (dx ** 2)
alpha_y = D0 * dt / (dy ** 2)

# Fonction M1(x, y, t) (exemple sinus pour illustrer)
def M1(x, y, t):
    return 0.1 * (1 + 0.5 * np.sin(2 * np.pi * x / Lx) * np.cos(2 * np.pi * y / Ly) * np.cos(0.01 * t))

# Initialisation
X, Y = np.meshgrid(np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny))
C = np.ones((Nx, Ny)) * C0
C_new = np.zeros((Nx, Ny))
normes_L2 = []

# Construction de la matrice tridiagonale (pour la diffusion implicite en 2D)
I = np.eye(Nx * Ny)
A = I * (1 + 2 * (alpha_x + alpha_y))

# Construction des décalages pour la matrice tridiagonale 2D
for i in range(Nx * Ny):
    if (i + 1) % Nx != 0:
        A[i, i + 1] = -alpha_x
        A[i + 1, i] = -alpha_x
    if i + Nx < Nx * Ny:
        A[i, i + Nx] = -alpha_y
        A[i + Nx, i] = -alpha_y

# Boucle temporelle
for n in range(Nt):
    t = n * dt
    M = M1(X, Y, t)
    B = C - k1 * M * C * dt  # Réaction explicite
    B_flat = B.flatten()
    C_new_flat = np.linalg.solve(A, B_flat)  # Diffusion implicite
    C_new = C_new_flat.reshape(Nx, Ny)
    C_new[0, :] = C_new[-1, :] = C_new[:, 0] = C_new[:, -1] = 0.0  # Conditions Dirichlet
    C[:] = C_new[:]
    normes_L2.append(np.sqrt(np.sum(C**2)))

# Affichage de la norme L2
plt.figure(figsize=(10, 6))
plt.plot(np.linspace(0, T, Nt), normes_L2, label="Norme L2")
plt.xlabel("Temps (s)")
plt.ylabel("Norme L2")
plt.title("Évolution de la norme L2 en 2D")
plt.grid(True)
plt.legend()
plt.show()

# Affichage de la concentration finale
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, C, levels=50, cmap="viridis")
plt.colorbar(label="Concentration")
plt.xlabel("Position x")
plt.ylabel("Position y")
plt.title("Concentration finale en 2D")
plt.show()


