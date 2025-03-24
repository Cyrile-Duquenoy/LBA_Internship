import numpy as np
import matplotlib.pyplot as plt

# Paramètres du problème
L = 1.0       # Longueur du domaine (carré)
Nx, Ny = 50, 50  # Nombre de mailles en x et y
Dx = L / Nx   # Pas spatial
Dy = L / Ny
Dt = 0.001    # Pas de temps
Tmax = 1.0    # Temps final
D = 0.1       # Coefficient de diffusion
r = 1.0       # Terme de réaction

# Création de la grille
x = np.linspace(0, L, Nx)
y = np.linspace(0, L, Ny)
X, Y = np.meshgrid(x, y)

# Initialisation de la concentration
U = np.exp(-((X-0.5)**2 + (Y-0.5)**2) / 0.01)  # Un pic gaussien au centre

# Schéma explicite en temps
def solve(U, Nx, Ny, Dx, Dy, Dt, D, r, Tmax):
    U_new = U.copy()
    t = 0.0
    while t < Tmax:
        U_old = U_new.copy()
        
        # Discrétisation en volumes finis avec conditions de Neumann
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                d2Udx2 = (U_old[i+1, j] - 2 * U_old[i, j] + U_old[i-1, j]) / Dx**2
                d2Udy2 = (U_old[i, j+1] - 2 * U_old[i, j] + U_old[i, j-1]) / Dy**2
                
                reaction = r * U_old[i, j] * (1 - U_old[i, j])  # Exemple : réaction logistique
                
                U_new[i, j] = U_old[i, j] + Dt * (D * (d2Udx2 + d2Udy2) + reaction)
        
        # Conditions de Neumann (flux nul sur les bords)
        U_new[0, :] = U_new[1, :]
        U_new[-1, :] = U_new[-2, :]
        U_new[:, 0] = U_new[:, 1]
        U_new[:, -1] = U_new[:, -2]

        t += Dt
    
    return U_new

# Résolution
U_final = solve(U, Nx, Ny, Dx, Dy, Dt, D, r, Tmax)

# Affichage du résultat
plt.imshow(U_final, extent=[0, L, 0, L], origin='lower', cmap='inferno')
plt.colorbar(label="Concentration")
plt.title("Équation de réaction-diffusion en volumes finis")
plt.xlabel("x")
plt.ylabel("y")
plt.show()


