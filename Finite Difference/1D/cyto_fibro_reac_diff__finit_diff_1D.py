import numpy as np
import matplotlib.pyplot as plt

N = 10
dx = 1/(N-1)

T = 10
dt = 0.01
t = np.arange(0, T + dt, dt)

D = 0.1
D_F = 3.78e-4
D_C = 5.62e-2

def is_cfl(lbda):
    if lbda <= 1/4:
        print("\n Condition CFL : ", True)
    else : 
        print("\n Condition CFL : ", False)
    return lbda <= 1/4

def lbda(D,dt,dx):
    return D * dt / (dx)**2

lbda_C = lbda(D_C, dt, dx)
lbda_F = lbda(D_F, dt, dx)

def cyto(x):
    return x*(1-x)

x = np.linspace(0,1,N)

def mat(N,lbda):
    d = np.ones(N)*(1 - 2*lbda)
    d[0] = (1 - lbda)
    d[-1] = d[0]
    d1 = np.ones(N-1)*(lbda)
    return np.diag(d) + np.diag(d1,1) + np.diag(d1,-1)

def cyto_init(x):
    U = []
    for i in range(len(x)):
        U.append(cyto(x[i]))
    return U

C = cyto_init(x)

plt.plot(x,C)
plt.title("Cyto. Initial")
plt.show()

def fibro(x):
    return 1-x
    
def fibro_init(x):
    U = []
    for i in range(len(x)):
        U.append(fibro(x[i]))
    return U
    
F = fibro_init(x)

plt.plot(x,F)
plt.title("Fibro. Initial")
plt.show()

def cyto_reac(C, dt, cyto_prod, cyto_death):
    C += dt*C*(cyto_prod - cyto_death)
    return C

def laplacian(C, dx):
    U = np.zeros_like(C)
    for i in range(1, len(C)-1):
        U[i] = (U[i-1] - 2*U[i] + U[i+1]) / (dx**2)
        U[0] = U[1]
        U[-1] = U[-2]
    return U

def fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi):
    return F + dt*(fibro_prod*F)*C  - chi*laplacian(C,dx) - fibro_death*F


# Fonction pour l'affichage côte à côte
def plot_F_C(F, C, i, t):
    plt.clf()  # Effacer la figure précédente
    plt.subplot(1, 2, 1)
    plt.plot(F, label="Fibroblastes", color='blue')
    plt.title(f"Fibroblastes (itération {i}/{len(t)})")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(C, label="Cytokines", color='red')
    plt.title(f"Cytokines (itération {i}/{len(t)})")
    plt.legend()

    plt.pause(0.5)  # Pause pour mise à jour dynamique

# PARAM.
cyto_prod = 1.2e-1
cyto_death = 1.2e-3
fibro_prod = 5.6
fibro_death = 4.8e-3
sat_f = 5.4e-1
chi = 0.003

if is_cfl(lbda_C) and is_cfl(lbda_F):
    A_C = mat(N, lbda_C)
    A_F = mat(N, lbda_F)
    plt.figure(figsize=(10, 5))
    for i in range(1,len(t)):
        C = np.dot(A_C, C)
        C = cyto_reac(C, dt, cyto_prod, cyto_death)
    
        F = np.dot(A_F, F)
        F = fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi)
    
        plot_F_C(F, C, i, t)  # Mise à jour du plot
    
    plt.show()  # Affichage final
    
    

    
