import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
from matplotlib import cm
import time

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
    x0, y0, z0 = N // 2, N // 4, N // 4
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

# Fonction pour créer des masques de voxels avec couleurs basées sur les valeurs absolues
def create_value_based_voxels(data, min_val=0):
    # S'assurer que les valeurs négatives sont exclues
    data_filtered = np.copy(data)
    data_filtered[data_filtered < min_val] = min_val
    
    # Définir les seuils à partir des valeurs réelles (non normalisées)
    # On utilise un seuil minimal (min_val) sans définir de seuil maximal
    # Pour éviter les valeurs nulles ou négatives
    data_positive = data_filtered[data_filtered > min_val]
    
    if len(data_positive) == 0:
        return [], [], 0, 0  # Pas de données positives
    
    # Définir les niveaux en fonction des valeurs réelles
    levels = 5  # Nombre de paliers de couleur
    
    min_data = np.min(data_positive) if len(data_positive) > 0 else min_val
    max_data = np.max(data)
    
    # Créer des seuils linéaires entre min_data et max_data
    if min_data == max_data:  # Si toutes les valeurs sont identiques
        thresholds = [min_data] * levels
    else:
        step = (max_data - min_data) / levels
        thresholds = [min_data + i * step for i in range(levels)]
    
    # Créer les masques pour chaque seuil
    masks = []
    for threshold in thresholds:
        mask = data > threshold
        masks.append(mask)
    
    # Créer un tableau de couleurs en fonction du niveau
    color_intensities = [0.2 + 0.8 * (i / (levels - 1)) for i in range(levels)]  # De 0.2 à 1.0
    
    return masks, color_intensities, min_data, max_data

# Fonction pour afficher les voxels colorés avec des valeurs absolues
def plot_combined_visualization(C, F, Coll, iteration, day):
    # Créer une figure avec deux rangées pour voxels et coupes
    fig = plt.figure(figsize=(18, 12))
    
    # Préparer les données pour les trois variables
    C_masks, C_colors, C_min, C_max = create_value_based_voxels(C)
    F_masks, F_colors, F_min, F_max = create_value_based_voxels(F)
    Coll_masks, Coll_colors, Coll_min, Coll_max = create_value_based_voxels(Coll)
    
    # Couleurs pour les voxels et les cartes de chaleur
    c_cmap = cm.Reds
    f_cmap = cm.Greens
    coll_cmap = cm.Blues
    
    # Première rangée: Visualisation des voxels
    # Sous-plot pour C (cytokines)
    ax1 = fig.add_subplot(2, 3, 1, projection='3d')
    for i, (mask, intensity) in enumerate(zip(C_masks, C_colors)):
        if np.any(mask):  # Vérifier si le masque contient des valeurs True
            ax1.voxels(mask, facecolors=colors.to_rgba(f'red', alpha=0.7*intensity), 
                      edgecolors=colors.to_rgba('darkred', alpha=0.1*intensity))
    
    ax1.set_title(f'Voxels - Cytokines \nMin: {C_min:.3e}, Max: {C_max:.3e}')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    
    # Sous-plot pour F (fibroblastes)
    ax2 = fig.add_subplot(2, 3, 2, projection='3d')
    for i, (mask, intensity) in enumerate(zip(F_masks, F_colors)):
        if np.any(mask):
            ax2.voxels(mask, facecolors=colors.to_rgba(f'green', alpha=0.7*intensity), 
                      edgecolors=colors.to_rgba('darkgreen', alpha=0.1*intensity))
    
    ax2.set_title(f'Voxels - Fibroblastes \nMin: {F_min:.3e}, Max: {F_max:.3e}')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    
    # Sous-plot pour Coll (collagène)
    ax3 = fig.add_subplot(2, 3, 3, projection='3d')
    for i, (mask, intensity) in enumerate(zip(Coll_masks, Coll_colors)):
        if np.any(mask):
            ax3.voxels(mask, facecolors=colors.to_rgba(f'blue', alpha=0.7*intensity), 
                      edgecolors=colors.to_rgba('darkblue', alpha=0.1*intensity))
    
    ax3.set_title(f'Voxels - Collagène \nMin: {Coll_min:.3e}, Max: {Coll_max:.3e}')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    
    # Deuxième rangée: Visualisation des coupes
    # Calculer les indices du centre de la grille
    center_idx = N // 2
    
    # Créer une grille 2D pour les plans XY, XZ et YZ
    x_2d, y_2d = np.meshgrid(x, y)
    x_2d_xz, z_2d_xz = np.meshgrid(x, z)
    y_2d_yz, z_2d_yz = np.meshgrid(y, z)
    
    # Sous-plot pour les coupes de C (cytokines)
    ax4 = fig.add_subplot(2, 3, 4, projection='3d')
    
    # Coupes selon les trois plans centrés
    xy_slice = C[:, :, center_idx]
    xz_slice = C[:, center_idx, :]
    yz_slice = C[center_idx, :, :]
    
    # Dessiner les trois plans
    surf1 = ax4.plot_surface(x_2d, y_2d, np.ones_like(x_2d) * z[center_idx], 
                            facecolors=c_cmap(xy_slice/np.max(C) if np.max(C) > 0 else xy_slice),
                            alpha=0.8, shade=False)
    surf2 = ax4.plot_surface(x_2d_xz, np.ones_like(x_2d_xz) * y[center_idx], z_2d_xz, 
                            facecolors=c_cmap(xz_slice/np.max(C) if np.max(C) > 0 else xz_slice),
                            alpha=0.8, shade=False)
    surf3 = ax4.plot_surface(np.ones_like(y_2d_yz) * x[center_idx], y_2d_yz, z_2d_yz, 
                            facecolors=c_cmap(yz_slice/np.max(C) if np.max(C) > 0 else yz_slice),
                            alpha=0.8, shade=False)
    
    ax4.set_title(f'Coupes - Cytokines')
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')
    sm = plt.cm.ScalarMappable(cmap=c_cmap, norm=plt.Normalize(C_min, C_max))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax4, shrink=0.5)
    
    # Sous-plot pour les coupes de F (fibroblastes)
    ax5 = fig.add_subplot(2, 3, 5, projection='3d')
    
    # Coupes selon les trois plans centrés
    xy_slice = F[:, :, center_idx]
    xz_slice = F[:, center_idx, :]
    yz_slice = F[center_idx, :, :]
    
    # Dessiner les trois plans
    surf1 = ax5.plot_surface(x_2d, y_2d, np.ones_like(x_2d) * z[center_idx], 
                            facecolors=f_cmap(xy_slice/np.max(F) if np.max(F) > 0 else xy_slice),
                            alpha=0.8, shade=False)
    surf2 = ax5.plot_surface(x_2d_xz, np.ones_like(x_2d_xz) * y[center_idx], z_2d_xz, 
                            facecolors=f_cmap(xz_slice/np.max(F) if np.max(F) > 0 else xz_slice),
                            alpha=0.8, shade=False)
    surf3 = ax5.plot_surface(np.ones_like(y_2d_yz) * x[center_idx], y_2d_yz, z_2d_yz, 
                            facecolors=f_cmap(yz_slice/np.max(F) if np.max(F) > 0 else yz_slice),
                            alpha=0.8, shade=False)
    
    ax5.set_title(f'Coupes - Fibroblastes')
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_zlabel('Z')
    sm = plt.cm.ScalarMappable(cmap=f_cmap, norm=plt.Normalize(F_min, F_max))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax5, shrink=0.5)
    
    # Sous-plot pour les coupes de Coll (collagène)
    ax6 = fig.add_subplot(2, 3, 6, projection='3d')
    
    # Coupes selon les trois plans centrés
    xy_slice = Coll[:, :, center_idx]
    xz_slice = Coll[:, center_idx, :]
    yz_slice = Coll[center_idx, :, :]
    
    # Dessiner les trois plans
    surf1 = ax6.plot_surface(x_2d, y_2d, np.ones_like(x_2d) * z[center_idx], 
                            facecolors=coll_cmap(xy_slice/np.max(Coll) if np.max(Coll) > 0 else xy_slice),
                            alpha=0.8, shade=False)
    surf2 = ax6.plot_surface(x_2d_xz, np.ones_like(x_2d_xz) * y[center_idx], z_2d_xz, 
                            facecolors=coll_cmap(xz_slice/np.max(Coll) if np.max(Coll) > 0 else xz_slice),
                            alpha=0.8, shade=False)
    surf3 = ax6.plot_surface(np.ones_like(y_2d_yz) * x[center_idx], y_2d_yz, z_2d_yz, 
                            facecolors=coll_cmap(yz_slice/np.max(Coll) if np.max(Coll) > 0 else yz_slice),
                            alpha=0.8, shade=False)
    
    ax6.set_title(f'Coupes - Collagène')
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    ax6.set_zlabel('Z')
    sm = plt.cm.ScalarMappable(cmap=coll_cmap, norm=plt.Normalize(Coll_min, Coll_max))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax6, shrink=0.5)
    
    plt.suptitle(f'Itération {iteration}, Jour {day:.2f}', fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.show()
    plt.pause(0.01)  # Pause pour permettre la mise à jour du graphique

# Initialisation
C = cyto_init(x, y, z)
F = fib_init(x, y, z)
Coll = coll_init()

# Simulation
if is_cfl(lbda_C) and is_cfl(lbda_F) and is_cfl(lbda_Coll):
    F_norm = [np.max(F)]
    C_norm = [np.max(C)]
    Coll_norm = [np.max(Coll)]
    
    print("Début de la simulation...")
    
    # Afficher les conditions initiales
    plot_combined_visualization(C, F, Coll, 0, 0)
    
    # Fréquence d'affichage (tous les X pas de temps)
    display_freq = 1000  # Ajustez cette valeur selon vos besoins
    
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
        
        # Afficher les résultats tous les display_freq pas de temps
        if i % display_freq == 0:
            print(f"Itération {i}/{len(t)}, Jour {day:.2f}")
            plot_combined_visualization(C, F, Coll, i, day)
        
        # Vérifier les valeurs négatives
        if is_negative(C) or is_negative(F) or is_negative(Coll):
            print(f"Negative Value : {i}, jour {day}")
            print(f"Minimum - C: {np.min(C)}, F: {np.min(F)}, Coll: {np.min(Coll)}")
    
    print("Simulation terminée.")
    
    # Afficher les résultats finaux
    plot_combined_visualization(C, F, Coll, len(t)-1, T)
    
    # Graphique de l'évolution des concentrations maximales
    plt.figure(figsize=(10, 6))
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

