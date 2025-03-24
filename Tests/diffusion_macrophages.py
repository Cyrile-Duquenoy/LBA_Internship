import numpy as np
import matplotlib.pyplot as plt

# Paramètres
L = 10.0        # Longueur de la région spatiale
T = 100.0         # Temps total de simulation
Nx = 100        # Nombre de points spatiaux
Nt = 500        # Nombre de pas de temps
k1 = 0.1        # Constante k1
M1 = 0.1        # Constante M1
D0 = 8.64*10**(-7)     # Coefficient de diffusion
C0 = 1.0        # Condition initiale de la concentration
dx = L / (Nx - 1) # Pas d'espace
dt = T / Nt # Pas de temps
alpha = D0 * dt / (dx ** 2)

# Initialisation des conditions initiales
C = np.ones(Nx) * C0  # Concentration initiale
C_new = np.zeros(Nx)  # Nouveau tableau pour la concentration mise à jour

# Liste pour stocker les concentrations à chaque pas de temps pour visualisation
concentrations = []

normes_L2 = []

# Boucle temporelle
for n in range(Nt):
    # Sauvegarde de la concentration à chaque étape
    concentrations.append(C.copy())

    # Mise à jour de la concentration à chaque point spatial
    for i in range(1, Nx - 1):
        C_new[i] = C[i] - k1 * M1 * C[i] * dt + D0 * (C[i+1] - 2 * C[i] + C[i-1]) * alpha
    
    # Appliquer les conditions aux limites (par exemple, C(0) et C(L) sont fixes)
    C_new[0] = C_new[-1] = 0.0  # Par exemple, condition de Dirichlet

    # Mise à jour de C
    C[:] = C_new[:]
    
    norme_L2 = np.sqrt(np.sum(C**2))
    normes_L2.append(norme_L2)

# Affichage des résultats
X = np.linspace(0, L, Nx)
plt.figure(figsize=(8, 6))
for i in range(0, Nt, Nt // 5):  # Affichage de 5 états intermédiaires
    plt.plot(X, concentrations[i], label=f'Time = {i * dt:.2f} s')
    
plt.xlabel('Position (x)')
plt.ylabel('Concentration (C)')
plt.title('Evolution de la concentration en fonction du temps')
plt.legend()
plt.grid(True)
plt.show()

# Plot de la norme L2
plt.figure(figsize=(8, 6))
plt.plot(np.linspace(0, T, Nt), normes_L2, label="Norme L2")
plt.xlabel('Temps (s)')
plt.ylabel('Norme L2')
plt.title('Evolution de la norme L2 de la concentration')
plt.grid(True)
plt.legend()
plt.show()






