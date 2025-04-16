import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp


# Paramètres
N = 10  # Taille de la grille 3D
L = 1.0
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
z = np.linspace(0, L, N)
X, Y, Z = np.meshgrid(x, y, z)
dx = L / (N - 1)
T = 20
dt = 0.001
t = np.arange(0, T + dt, dt)

# Coefficients de diffusion
D_C = 2.58e-2
D_F = 1.46e-7
D_Coll = 4.59e-8
rho = 1e-4

def lbda(D, dt, dx): return D * dt / (dx**2)

lbda_C = lbda(D_C, dt, dx)
lbda_F = lbda(D_F, dt, dx)
lbda_Coll = lbda(D_Coll, dt, dx)

def is_cfl(lbda): return lbda <= 1/6  # Condition CFL pour la 3D
def is_negative(U): return np.min(U) < 0

# Paramètres de réaction
fib_val = 1e-2
fib_prod = 1
fib_death = 0.1
cyto_prod = 1
cyto_death = 1
coll_prod = 0.6
coll_death = 0.6
chi = 2
sat_F = 1
coll_sat = 1
chi_coll = 1e-10

def cyto(x, y, z):
    return (x * (1 - x) + y * (1 - y) + z * (1 - z)) * 1.13e-1

def cyto_init(x, y, z):
    U = np.zeros((N, N, N))
    x0, y0, z0 = N // 2, N // 2, N // 2
    rad = N // 2
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if (i - x0)**2 + (j - y0)**2 + (k - z0)**2 < rad**2:
                    U[i, j, k] = cyto(x[i], y[j], z[k])
    return U

def fib(x, y, z):
    return fib_val * x * (1 - y) * (0.5 + 0.5 * z)

def fib_init(x, y, z):
    U = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                U[i, j, k] = fib(x[i], y[j], z[k])
    return U

def coll_init():
    return np.zeros((N, N, N))

def laplacian_3d(U, dx):
    L = np.zeros_like(U)
    
    # Intérieur du domaine
    L[1:-1, 1:-1, 1:-1] = (
        U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
        U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
        U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
        6 * U[1:-1, 1:-1, 1:-1]
    ) / (dx**2)
    
    # Neumann Conditions
    # Face x=0
    L[0, :, :] = L[1, :, :]
    # Face x=L
    L[-1, :, :] = L[-2, :, :]
    # Face y=0
    L[:, 0, :] = L[:, 1, :]
    # Face y=L
    L[:, -1, :] = L[:, -2, :]
    # Face z=0
    L[:, :, 0] = L[:, :, 1]
    # Face z=L
    L[:, :, -1] = L[:, :, -2]
    
    return L

# Réactions
def cyto_reac(C, F, cyto_prod, cyto_death):
    return cyto_prod * F * C - cyto_death * C

def fib_reac(F, C, fib_prod, fib_death, chi, lapC):
    term = fib_prod * C * (1 - F / sat_F)
    death = fib_death * F
    chemotaxis = chi * F * lapC
    return term - death - chemotaxis

def coll_reac(Coll, F, coll_prod, coll_death, coll_sat, lapF):
    term = coll_prod * F * (1 - Coll / coll_sat)
    death = coll_death * Coll
    chemotaxis = chi_coll * Coll * lapF
    return term - death - chemotaxis

# Fonction pour les schémas de diffusion explicite 3D
def diffuse_explicit_3d(U, lbda):
    V = U.copy()
    V[1:-1, 1:-1, 1:-1] = U[1:-1, 1:-1, 1:-1] + lbda * (
        U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
        U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
        U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
        6 * U[1:-1, 1:-1, 1:-1]
    )
    
    # Neumann Conditions
    V[0, :, :] = V[1, :, :]
    V[-1, :, :] = V[-2, :, :]
    V[:, 0, :] = V[:, 1, :]
    V[:, -1, :] = V[:, -2, :]
    V[:, :, 0] = V[:, :, 1]
    V[:, :, -1] = V[:, :, -2]
    
    return V

# Initialisation
C = cyto_init(x, y, z)
F = fib_init(x, y, z)
Coll = coll_init()

# Simulation
if is_cfl(lbda_C) and is_cfl(lbda_F) and is_cfl(lbda_Coll):
    F_norm = [np.max(F)]
    C_norm = [np.max(C)]
    Coll_norm = [np.max(Coll)]

    
    # Collecter les résultats pour les animations
    C_results = [C.copy()]
    F_results = [F.copy()]
    Coll_results = [Coll.copy()]
    
    print("Début de la simulation...")
    
    for i in range(1, len(t)):
        day = i * dt
        
        # Diffusion et réaction pour C (cytokines)
        C = diffuse_explicit_3d(C, lbda_C)
        lapC = laplacian_3d(C, dx)
        C += dt * cyto_reac(C, F, cyto_prod, cyto_death)
        C_norm.append(np.max(C))
        
        # Diffusion et réaction pour F (fibroblastes)
        F = diffuse_explicit_3d(F, lbda_F)
        lapF = laplacian_3d(F, dx)
        F += dt * fib_reac(F, C, fib_prod, fib_death, chi, lapC)
        F_norm.append(np.max(F))
        
        # Diffusion et réaction pour Coll (collagène)
        Coll = diffuse_explicit_3d(Coll, lbda_Coll)
        Coll += dt * coll_reac(Coll, F, coll_prod, coll_death, coll_sat, lapF)
        Coll_norm.append(np.max(Coll))
        
        # Vérifier les valeurs négatives
        if is_negative(C) or is_negative(F) or is_negative(Coll):
            print(f"Negative Value : {i}, jour {day}")
            print(f"Minimum - C: {np.min(C)}, F: {np.min(F)}, Coll: {np.min(Coll)}")
            
    
    print("Simulation terminée.")
    
    
    # Résultats
    plt.plot(t, F_norm, label='Fibro.')
    plt.plot(t, C_norm, label='Cyto.')
    plt.plot(t, Coll_norm, label='Coll.')
    plt.xlabel("Temps (j)")
    plt.ylabel("Valeur max.")
    plt.title("Évolution des concentrations maximales")
    plt.legend()
    plt.grid()
    plt.show()

    
else:
    print("Erreur: La condition CFL n'est pas satisfaite.")
    print(f"lambda_C = {lbda_C}, lambda_F = {lbda_F}, lambda_Coll = {lbda_Coll}")
    print("Ces valeurs doivent être inférieures à 1/6 pour garantir la stabilité en 3D.")

