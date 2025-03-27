import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

N = 10
dx = 1/(N-1)
L = 1
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y) 

T = 10
dt = 0.1
t = np.arange(0, T + dt, dt)

D_F = 1.46e-4
D_C = 2.58e-2

lbda_C =  D_C * dt /(dx**2)
lbda_F =  D_F * dt /(dx**2)

rho = 1e-4

# %%

def is_cfl(lbda:float):
    if lbda <= 1/4:
        print("\n Condition CFL : ", True)
    else : 
        print("\n Condition CFL : ", False)
    return lbda <= 1/4

def lbda(D:float,dt:float,dx:float)->float:
    return D * dt / (dx)**2

# %%
def diffusion_2_neumann_2d(U:list, lbda:float)->list:
    N = len(U)
    U_new = np.zeros_like(U)
    for i in range(1,len(U)-1):
        for j in range(1,len(U[i])-1):
            U_new[i,j] = (1 - 4*lbda)*U[i, j] + lbda*( U[i+1, j] + U[i-1, j] + U[i, j+1] + U[i, j-1] )
    
    # Gestion des bords
    U_new[0,:] = U[1,:] 
    U_new[-1,:] = U[-2,:]
    U_new[:,0] = U[:,1]
    U_new[:,-1] = U[:,-2]

    # Gestion des coins
    U[0,0] = U[1,1]
    U[N-1,N-1] = U[N-2,N-2]
    U[0,N-1] = U[1,N-2]
    U[N-1,0] = U[N-2,1]
    
    return U_new

# %%

def cyto(x:float, y:float)->float:
    return x*(1-y)*1.13e-1

def cytok_init(x:list,y:list,rho:float):
    x0, y0 = len(x) // 4, len(y) // 4
    rad = min(len(x),len(y)) // 4
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            U[i,j] = cyto(x[i], y[j])
    return U

C = cytok_init(x,y,rho)

# %%



def fibro(x:float, y:float)->float:
    return 1e-4

def fib_init(x:list, y:list)->list:
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            U[i,j] = fibro(x[i], y[j])
    return U


F = fib_init(x, y)

# %%

def chemotaxis_flux(F, C, chi, dx, N):
    """
    Calcule ∇⋅(χ F ∇C) en différences finies avec condition de Neumann homogène.

    Retourne :
    - Matrice (N, N) représentant ∇⋅( χF ∇C)
    """
    chi_F = chi * F
    div_chiF_gradC = np.zeros((N, N))
    for i in range(1, N-1):
        for j in range(1, N-1):
            dCdx = (C[i+1, j] - C[i-1, j]) / (2 * dx)
            dCdy = (C[i, j+1] - C[i, j-1]) / (2 * dx)

            d_chiF_dCdx = (chi_F[i+1, j] * dCdx - chi_F[i-1, j] * dCdx) / (2 * dx)
            d_chiF_dCdy = (chi_F[i, j+1] * dCdy - chi_F[i, j-1] * dCdy) / (2 * dx)

            div_chiF_gradC[i, j] = d_chiF_dCdx + d_chiF_dCdy
            
    # Neumann conditions
    div_chiF_gradC[0, :] = div_chiF_gradC[1, :]
    div_chiF_gradC[-1, :] = div_chiF_gradC[-2, :]
    div_chiF_gradC[:, 0] = div_chiF_gradC[:, 1]
    div_chiF_gradC[:, -1] = div_chiF_gradC[:, -2]
    
    return div_chiF_gradC

# %% Param
fibro_prod = 1
fibro_death = 1
cyto_prod = 1
cyto_death = 1
sat_f =1
chi = 1

# %% Reaction Terms

def fibro_reaction(F, C, fibro_prod, chi, fibro_death, sat_f, N):
    return fibro_prod*C*(1 - F/sat_f) - diffusion_2_neumann_2d(C, lbda_C)  - fibro_death*F

def cyto_reaction(C, F, cyto_prod, cyto_death):
    return cyto_prod*C - cyto_death*C

# %%

min_val = 0
max_val_C = C.max()
max_val_F = 3*(F).max()  # Mets la valeur maximale souhaitée

fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 ligne, 2 colonnes

im1 = axes[0].imshow(C, cmap='viridis', origin='lower', extent=[0, L, 0, L], vmin=min_val, vmax=max_val_C)
axes[0].set_title("Cytokines initiales")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
fig.colorbar(im1, ax=axes[0], label="Concentration de cytokines")

im2 = axes[1].imshow(F, cmap='viridis', origin='lower', extent=[0, L, 0, L], vmin=min_val, vmax=max_val_F)
axes[1].set_title("Fibroblastes Initiales")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
fig.colorbar(im2, ax=axes[1], label="Concentration de fibroblastes")

plt.tight_layout()  # Ajuste l'espacement pour éviter le chevauchement
plt.show()

# %%

cytok_norm = []
fib_norm = []

if is_cfl(lbda_C) and is_cfl(lbda_F):
    for n in range(len(t)):
        
        C = diffusion_2_neumann_2d(C, lbda_C)
        C += cyto_reaction(C, F, cyto_prod, cyto_death)
        F = diffusion_2_neumann_2d(F, lbda_F)
        F += dt*fibro_reaction(F, C, fibro_prod, chi, fibro_death, sat_f, N)
        
        cytok_norm.append(C.max())
        fib_norm.append(F.max())

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 ligne, 2 colonnes
        im1 = axes[0].imshow(C, cmap='viridis', origin='lower', extent=[0, L, 0, L], vmin=min_val, vmax=max_val_C)
        axes[0].set_title(f"Cytokines à t = {t[n]:.2f} j")
        axes[0].set_xlabel("x")
        axes[0].set_ylabel("y")
        fig.colorbar(im1, ax=axes[0], label="Concentration de cytokines")
        
        im2 = axes[1].imshow(F, cmap='viridis', origin='lower', extent=[0, L, 0, L], vmin=min_val, vmax=max_val_F)
        axes[1].set_title(f"Fibroblastes à t = {t[n]:.2f} j")
        axes[1].set_xlabel("x")
        axes[1].set_ylabel("y")
        fig.colorbar(im2, ax=axes[1], label="Concentration de fibroblastes")
        
        plt.tight_layout()  # Ajuste l'espacement pour éviter le chevauchement
        
        plt.show()
plt.show()
        
# %%

plt.figure()
plt.plot(t[:len(cytok_norm)], cytok_norm, label="\u03C3 (Cytokine)/Jour")
plt.xlabel("Temps (j)")
plt.ylabel("Cytokine moyenne / Unité")
#plt.yscale('log')
#plt.xscale('log')
plt.title("Cytokine")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(t[:len(fib_norm)], fib_norm, label="\u03C3 (Fibroblaste)/Jour")
plt.xlabel("Temps (j)")
plt.ylabel("Fibroblaste moyenne / Unité")
#plt.yscale('log')
#plt.xscale('log')
plt.title("Fibroblaste")
plt.legend()
plt.grid()
plt.show()




