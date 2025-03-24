import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg
import matplotlib.pyplot as plt

'''######################################################
   #######  Simple diffusion 2D. ########################
   #######  Ce code utilise des sparses matrix. #########
   #######  Prochaine étape : ###########################
   #######      ==> Inclure le terme de réaction ########
   ######################################################
'''

# Paramètres
N = 50
M = 200
L = 1
T = 200
dx = L / (N - 1)
dt = 0.01

D = 1.47*1e-6  # Coeff. de diffusion
lbda = D * dt / dx**2

rho = 1e5  # Densité initiale des fibroblastes
rad = N // 2
x0, y0 = N // 2, N // 2

# Initialisation de U
# Pour Concentrer les fibroblastes au temp t=0 dans périmètre défini
def u_init():
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad: 
                U[i, j] = rho
    return U

U = u_init().flatten()

t = np.arange(0, T + dt, dt)

# Construction des sparses matrix
def Mat_sparse(N):
    size = N * N
    diag_main = np.ones(size) * (1 - 4 * lbda)
    diag_x = np.ones(size - 1) * lbda
    diag_y = np.ones(size - N) * lbda
    
    for i in range(1, N):
        diag_x[i * N - 1] = 0  # Gestion des bords
    
    A = sp.diags([diag_main, diag_x, diag_x, diag_y, diag_y], [0, -1, 1, -N, N], format='csr')
    return A

A = Mat_sparse(N)
l2_norms = []

plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label="Concentration (cellules/cm²)")
plt.title("Distribution initiale des fibroblastes hépatiques")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# Évolution temporelle
for n in range(len(t)):
    U = A @ U  # Produit matrice creuse/vecteur
    l2_norms.append(np.linalg.norm(U, ord=2))
    
    if n % 4000 == 0:
        plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
        plt.colorbar(label="Concentration (cellules/cm²)")
        plt.title(f"Évolution à t = {t[n]:.2f} j")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.pause(0.5)

plt.show()

# Plot de la norme L2
plt.figure()
plt.plot(t[:len(l2_norms)], l2_norms, label="Norme L2 de la concentration")
plt.xlabel("Temps (j)")
plt.ylabel("Norme L2")
plt.yscale('log')
#plt.xscale('log')
plt.title("Évolution de la norme L2 de la concentration")
plt.legend()
plt.grid()
plt.show()


