import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# Paramètres
N = 50
L = 1
T = 200
dx = L / (N - 1)
dt = 0.0001

D = 1.47e-6  # Coeff. de diffusion
lbda = D * dt / dx**2

rho = 1.2e-2
rad = N // 2
x0, y0 = N // 2, N // 2

lbda_Ef = 2.5*1e-1 
lbda_fHA = 2.5*1e-3
lbda_fE = 5*1e-4 
lbda_mfT = 0.12
lbda_mfG = 0.12
E0 = 0.1
H_A = 1*1e-4
K_HA = 2*1e-3
T_b = 1.52*1e-9
K_Tb = 1*1e-10
G = 3.07*1e-10
K_G = 1.5*1e-8
d_m = 1.66*1e-2
I13 = 1.13*1e-9
K_I13 = 2*1e-7
E = 1*1e-6
K_E = 0.1
d_f = 1.66*1e-2

#################################################

def is_cfl(dt,dx):
    return D*dt <= (dx**2)/4

print(is_cfl(dt,dx))

###############################################################################

def u_init():
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad:
                U[i, j] = rho
    return U

U = u_init().flatten()
t = np.arange(0, T + dt, dt)

###############################################################################

def Mat_sparse(N):
    size = N * N
    diag_main = np.ones(size) * (1 - 4 * lbda)
    diag_x = np.ones(size - 1) * lbda
    diag_y = np.ones(size - N) * lbda
    
    for i in range(1, N):
        diag_x[i * N - 1] = 0  # Gestion des bords
    
    A = sp.diags([diag_main, diag_x, diag_x, diag_y, diag_y], [0, -1, 1, -N, N], format='csr')

    # Conditions de Neumann
    for i in range(N):  
        for j in range(N):
            index = i * N + j  # Index global
            if i == 0:
                A[index, index + N] += lbda  # Bord bas
            if i == N - 1:
                A[index, index - N] += lbda  # Bord haut
            if j == 0:
                A[index, index + 1] += lbda  # Bord gauche
            if j == N - 1:
                A[index, index - 1] += lbda  # Bord droit
    
    return A

A = Mat_sparse(N)

plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label="Concentration (cellules/cm²)")
plt.title("Distribution initiale des fibroblastes")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

############################# Reaction ########################################

def R(u):
    I1 = lbda_Ef * E0
    I2 = (H_A / (K_HA + H_A)) * u
    I3 = lbda_fE * (T_b / (K_Tb + T_b) + I13 / (K_I13 + I13)) * (E / (K_E + E)) * u
    I4 = (lbda_mfT * (T_b / (K_Tb + T_b)) + lbda_mfG * (G / (K_G + G))) * u
    I5 = d_f * u
    return I1 + I2 + I3 - I4 - I5

l2_norms = []
for n in range(1,len(t)):
    U = A @ U  # Diffusion
    U += dt * R(U)  # Réaction
    
    l2_norms.append(np.linalg.norm(U, ord=2))
    
    if n % 100000 == 0:
        plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
        plt.colorbar(label="Concentration (cellules/cm²)")
        plt.title(f"Évolution à t = {t[n]:.2f} j")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.pause(0.2)
plt.show()

plt.figure()
plt.plot(t[:len(l2_norms)], l2_norms, label="Norme L2 de la concentration")
plt.xlabel("Temps (j)")
plt.ylabel("Norme L2")
plt.yscale('log')
plt.title("Évolution de la norme L2 de la concentration")
plt.legend()
plt.grid()
plt.show()


