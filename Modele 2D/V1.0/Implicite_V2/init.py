import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

rho = 1.2e-2
D_C = 1.08e-2

def fib_init(N, x0, y0, rad, rho, sigma=2, scale=1.0, occupancy=0.06):
    # Générer une matrice avec des valeurs aléatoires uniformes
    random_matrix = np.random.rand(N, N)
    
    # Appliquer un filtre gaussien pour lisser la distribution
    fib_matrix = gaussian_filter(random_matrix, sigma=sigma)
    
    # Normaliser les valeurs entre 0 et 1
    fib_matrix = fib_matrix / np.max(fib_matrix)
    
    # Seuiller pour ne garder que les valeurs les plus élevées correspondant à 5% de la surface
    threshold = np.percentile(fib_matrix, 100 * (1 - occupancy))
    fib_matrix = np.where(fib_matrix >= threshold, 1.0, 0.0)
    
    # Mise à l'échelle
    fib_matrix *= scale
    '''
    for i in range(N):
        for j in range(N):
            fib_matrix[i,j] = 1e-6
    '''
    
    return fib_matrix


def cytok_init(N, sigma=2, scale=1.0, occupancy=1.):
    # Générer une matrice avec des valeurs aléatoires uniformes
    random_matrix = np.random.rand(N, N)
    
    # Appliquer un filtre gaussien pour lisser la distribution
    cytok_matrix = gaussian_filter(random_matrix, sigma=sigma)
    
    # Normaliser les valeurs entre 0 et 1
    cytok_matrix = cytok_matrix / np.max(cytok_matrix)
    
    # Mise à l'échelle
    cytok_matrix *= scale
    
    return cytok_matrix
'''



def cytok_init(N,rho):
    x0, y0 = N // 4, N // 4
    rad = N // 2
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad: 
                U[i, j] = rho
    return U



'''

def initialisation(N, rho, L, rad):
    rad=N//2
    C0=cytok_init(N,rho).flatten()
    F0=fib_init(N, N//2, N//2, rad, rho).flatten()
    
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 ligne, 2 colonnes

    im1 = axes[0].imshow(C0.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
    axes[0].set_title("Cytokines initiales en distribution Gaussienne")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im1, ax=axes[0], label="Concentration de cytokines")

    im2 = axes[1].imshow(F0.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
    axes[1].set_title("Fibroblastes initiaux en distribution Gaussienne")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    fig.colorbar(im2, ax=axes[1], label="Concentration de fibroblastes")

    plt.tight_layout()  # Ajuste l'espacement pour éviter le chevauchement
    plt.show()
    
    return C0,F0

if __name__ == "__main__":
    I=initialisation(50, 1.2 * 10e-2, 1)


