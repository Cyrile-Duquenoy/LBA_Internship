import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg
import matplotlib.pyplot as plt


'''######################################################
   #######  Simple diffusion 2D. ########################
   #######  Ce code utilise des matrices pleines. #######
   #######  Prochaine étape : ###########################
   #######      ==> Sparses Matrix ######################
   ######################################################
'''

# Paramètres
N = 50
M = 400
L = 1
T = 400
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)
dx = L / (N - 1)
dt = 0.01

def CFL_condition(dx:float,dt:float):
    return True

# Vecteur de temps défini
t = np.arange(0, T + dt, dt)

D = 1e-3  # Coeff. de diffusion
lbda = D * dt / dx**2

rho = 1e5  # Densité initiale des fibroblastes
rad = N // 2
x0, y0 = N // 2, N // 2

# Initialisation de U
def u_init() -> np.ndarray:
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad:
                U[i, j] = rho
    return U

U = u_init()

plt.imshow(U, cmap='viridis', origin='lower', extent=[0, L, 0, L])
plt.colorbar(label="Concentration (cellules/cm²)")
plt.title("Distribution initiale des fibroblastes hépatiques")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

U = U.flatten()  # Aplatir la matrice en un vecteur


# Matrices de discrétisation
def Mat_x(N: int) -> np.ndarray:
    size = N * N
    d = np.ones(size) * (1 - 4 * lbda)
    d1 = np.ones(size - 1) * (-lbda)
    for i in range(1, N):
        d1[i * N - 1] = 0  # Correction pour éviter les erreurs aux bords
    return np.diag(d) + np.diag(d1, 1) + np.diag(d1, -1)

def Mat_y(N: int) -> np.ndarray:
    size = N * N
    d = np.ones(size - N) * lbda
    return np.diag(d, -N) + np.diag(d, N)

Ax = Mat_x(N)
Ay = Mat_y(N)
A = Ax + Ay

l2_norms = []  # Stocke la norme L2 de U à chaque instant

# Évolution temporelle
for n in range(len(t)):  # Parcourir tout le vecteur temps
    U = A @ U  # Mise à jour
    l2_norms.append(np.linalg.norm(U, ord=2))  # Norme L2 de la concentration

    if n % 4000 == 0:  # Afficher tous les 10 pas de temps
        plt.imshow(U.reshape(N, N), cmap='viridis', origin='lower', extent=[0, L, 0, L])
        plt.colorbar(label="Concentration (cellules/cm²)")
        plt.title(f"Évolution à t = {t[n]:.2f} j")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.pause(0.5)
plt.show()

##### Plot de la norme L2 de U #########
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


