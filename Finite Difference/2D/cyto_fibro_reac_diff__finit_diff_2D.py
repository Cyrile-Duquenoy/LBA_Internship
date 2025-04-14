import numpy as np
import matplotlib.pyplot as plt

N = 10
L = 1. 
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)
dx = L / (N - 1)
T = 20
dt = 0.1
t = np.arange(0, T + dt, dt)

D_C = 2.58e-2
D_F = 1.46e-7
D_Coll = 4.59e-8
rho = 1e-4

def lbda(D:float, dt:float, dx:float)->float:
    return D * dt / (dx**2)

lbda_C = lbda(D_C, dt, dx)
lbda_F = lbda(D_F, dt, dx)
lbda_Coll = lbda(D_Coll, dt, dx)

def is_cfl(lbda):
    return lbda <= 1/4

def is_negative(U:list)->bool:
    return np.min(U) < 0

#%%

fib = 1e-2
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

#%%
def cyto(x:float, y:float, N)->float:
    return 1.13e-1

def cyto_init(x:list,y:list,rho:float):
    x0, y0 = len(x) // 2, len(y) // 4
    rad = min(len(x),len(y)) // 2
    U = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (i - x0) ** 2 + (j - y0) ** 2 < rad:
                U[i,j] = cyto(x[i], y[j], N)
    return U

C = cyto_init(x, y, rho)

def cyto_reac(C, F, cyto_prod, cyto_death):
    return cyto_prod*F*C - cyto_death*C

#%%

def laplacian(U:list, dx:float)->list:
    U_new = np.zeros_like(U)
    for i in range(1, len(U)-1):
        for j in range(1, len(U[i])-1):
            U_new[i, j] = (U[i-1, j] + U[i+1, j] - 4*U[i, j] + U[i, j+1] + U[i, j-1]) / (dx**2)
    # Gestion des bords
    U_new[0,:] = U[1,:] 
    U_new[-1,:] = U[-2,:]
    U_new[:,0] = U[:,1]
    U_new[:,-1] = U[:,-2]

    # Gestion des coins
    U_new[0,0] = U_new[1,1]
    U_new[N-1,N-1] = U_new[N-2,N-2]
    U_new[0,N-1] = U_new[1,N-2]
    U_new[N-1,0] = U_new[N-2,1]

    return U_new

#%%
def fib(x:float, y:float, fib=fib)->float:
    return fib*x*(1-y)

def fib_init(x:list, y:list):
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            U[i,j] = fib(x[i], y[j])
    return U

F = fib_init(x, y)

def fib_reac(F, C, fib_prod, fib_death, chi, N, dx):
    return fib_prod*C*(1 - F/sat_F) - fib_death*F - chi*F*laplacian(C, dx)

#%%

def coll(x:float, y:float)->float:
    return 0 

def coll_init(x:list, y:list)->list:
    U = np.zeros((N, N)) 
    for i in range(len(x)):
        for j in range(len(y)):
            U[i, j] = coll(x[i], y[j])
    return U

def coll_reac(Coll, F, coll_prod, coll_death, sat_coll, N, dx):
    return coll_prod*F*(1 - Coll/sat_coll) - coll_death*Coll - chi_coll*Coll*laplacian(F, dx)
    
Coll = coll_init(x, y)

#%%

def mat_diff(N:int, lbda:float)->list:
    size = N**2
    d = np.ones(size)*(1 - 4*lbda)
    d1 = np.ones(size - 1)*(lbda)
    d2 = np.ones(size - N)*(lbda)
    return np.diag(d) + np.diag(d1, 1) + np.diag(d1, -1) + np.diag(d2, N) + np.diag(d2, -N)

if(is_cfl(lbda_C) and is_cfl(lbda_F) and is_cfl(lbda_Coll)):
    L_C = mat_diff(N, lbda_C)
    L_F = mat_diff(N, lbda_F)
    L_Coll = mat_diff(N, lbda_Coll)
    
    F_norm = [np.max(F)]
    C_norm = [np.max(C)]
    Coll_norm = [np.max(Coll)]
    
    vmin_value = 0
    for i in range(1,len(t)):
        day = i*dt
        
        C = C.flatten()
        C = np.dot(L_C, C)
        C = C.reshape((N, N))
        C += dt*cyto_reac(C, F, cyto_prod, cyto_death)
        C_norm.append(np.max(C))
        
        F = F.flatten()
        F = np.dot(L_F, F)
        F = F.reshape((N, N))
        F += dt*fib_reac(F, C, fib_prod, fib_death, chi, N, dx)
        F_norm.append(np.max(F))
        
        Coll = Coll.flatten()
        Coll = np.dot(L_Coll, Coll)
        Coll = Coll.reshape((N, N))
        Coll += dt*coll_reac(Coll, F, coll_prod, coll_death, coll_sat, N, dx)
        Coll_norm.append(np.max(Coll))
        
        print("Iteration ", i)
        print("F has a negative value : ", is_negative(F))
        print("C has a negative value : ", is_negative(C))
        print("Coll has a negative value : ", is_negative(Coll))
        print("\n")
        
        # Création de la figure avec deux sous-graphes côte à côte
        fig = plt.figure(figsize=(10, 10))
            
        # --- Subplot 1 : C ---
        ax1 = fig.add_subplot(1, 3, 1, projection='3d')
        surf1 = ax1.plot_surface(X, Y, C, cmap='viridis', vmin=vmin_value)
        fig.colorbar(surf1, ax=ax1, shrink=0.3)
        ax1.set_title(f"Cytokine (C) - day {day}")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_zlabel("C")
        
        # --- Subplot 2 : C ---
        ax1 = fig.add_subplot(1, 3, 2, projection='3d')
        surf1 = ax1.plot_surface(X, Y, F, cmap='plasma', vmin=vmin_value)
        fig.colorbar(surf1, ax=ax1, shrink=0.3)
        ax1.set_title(f"Fibroblaste (F) - day {day}")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_zlabel("F")
        
        # --- Subplot 3 : Coll ---
        ax1 = fig.add_subplot(1, 3, 3, projection='3d')
        surf1 = ax1.plot_surface(X, Y, Coll, cmap='plasma', vmin=vmin_value)
        fig.colorbar(surf1, ax=ax1, shrink=0.3)
        ax1.set_title(f"Collagène (C) - day {day}")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_zlabel("Coll") 
        #plt.pause(0.1)
        plt.show()
        
plt.plot(t, F_norm, label='F_norm')
plt.plot(t, C_norm, label='C_norm')
plt.plot(t, Coll_norm, label='Coll_norm')
plt.show()
    
