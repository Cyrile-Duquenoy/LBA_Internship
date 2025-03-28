import numpy as np
import matplotlib.pyplot as plt

# %% PARAM.

cyto_prod = 1
cyto_death = 1
fibro_prod = 1
fibro_death = 0.1
sat_f = 1
chi = 2
coll_prod = 0.8
coll_death = 0.5
sat_coll = 1

cytok_norm = []
fib_norm = []
coll_norm = []
Index = []

#%%

# Space Discret.
N = 10
dx = 1/(N-1)
L = 1
x = np.linspace(0,L,N)

# Time Discret
T = 70
dt = 0.1
t = np.arange(0, T + dt, dt)

# Diffusion Ceoff.
D_F = 1.46e-7
D_C = 2.58e-2
D_Coll = 4.59e-13

def is_cfl(lbda):
    if lbda <= 1/2:
        print("\n Condition CFL : ", True)
    else : 
        print("\n Condition CFL : ", False)
    return lbda <= 1/2

def lbda(D,dt,dx):
    return D * dt / (dx)**2

lbda_C = lbda(D_C, dt, dx)
lbda_F = lbda(D_F, dt, dx)
lbda_Coll = lbda(D_Coll, dt, dx)

# %%
def mat(N,lbda):
    d = np.ones(N)*(1 - 2*lbda)
    d[0] = (1 - lbda)
    d[-1] = d[0]
    d1 = np.ones(N-1)*(lbda)
    return np.diag(d) + np.diag(d1,1) + np.diag(d1,-1)

#%%

def cyto(x):
    return x*(1-x)*1.13e-1

def cyto_init(x):
    U = []
    for i in range(len(x)):
        U.append(cyto(x[i]))
    return U

def cyto_reac(C, F, dt, cyto_prod, cyto_death):
    C += dt*C*F*cyto_prod - dt*C*cyto_death
    return C

C = cyto_init(x)

plt.plot(x,C)
plt.title("Cyto. Initial")
plt.show()

#%%

def fibro(x):
    return 1e-4
    
def fibro_init(x):
    U = []
    for i in range(len(x)):
        U.append(fibro(x[i]))
    return U

def fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi):
    return F + dt*( fibro_prod*C*(1 - F/sat_f) - chi*laplacian(C,dx) - fibro_death*F )
    
F = fibro_init(x)

plt.plot(x,F)
plt.title("Fibro. Initial")
plt.show()

#%%

def coll(x):
    return 0

def coll_init(x):
    U = []
    for i in range(len(x)):
        U.append(coll(x[i]))
    return U

def coll_reac(Coll, F, dt, coll_prod, coll_death):
    return Coll + dt*(coll_prod*F*(1 - Coll/sat_coll) - laplacian(F, dx) - coll_death*Coll)

Coll = coll_init(x)

# %%

def laplacian(U, dx):
    U_new = np.zeros_like(U)
    for i in range(1, len(C)-1):
        U_new[i] = (U_new[i-1] - 2*U_new[i] + U_new[i+1]) / (dx**2)
    U_new[0] = U_new[1]
    U_new[-1] = U_new[-2]
    return U

# %% Plot Function

def plot_F_C(F, C, Coll, i, t, max_val):
    plt.clf()  # Suppr. Previous Figure
    min_val = 0
    plt.subplot(1, 3, 1)
    plt.plot(F, label="Fibroblastes", color='blue')
    plt.plot(C, label="Cytokines", color='red')
    plt.plot(Coll, label="Collagène", color='green')
    plt.title(f"Fibroblastes, Cytokines & Collagene (itération {i}/{len(t)-1})")
    plt.ylim(min_val)
    plt.legend()
    plt.pause(0.01)

# %% Time Iterations

if is_cfl(lbda_C) and is_cfl(lbda_F): # Verify CFL conditions

    # Discret Matrix
    A_C = mat(N, lbda_C)
    A_F = mat(N, lbda_F)
    A_Coll = mat(N, lbda_Coll)
    
    # Append Norm
    cytok_norm.append(max(C))
    fib_norm.append(max(F))
    coll_norm.append(max(Coll))
    Index.append(np.argmax(F)/N)
    
    plt.figure(figsize=(10, 5))
    max_val = max(max(C), max(F))
    
    # Time Iteration
    for i in range(1,len(t)):
        C = np.dot(A_C, C)
        C = cyto_reac(C, F, dt, cyto_prod, cyto_death)
        '''
        # Repetitive Trauma Simulation
        if i%100 == 0:
            C+=cyto_init(x)
        '''
        m_c = min(C)
        
        cytok_norm.append(max(C))
        fib_norm.append(max(F))
    
        F = np.dot(A_F, F)
        F = fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi)
        m_f = min(F)
        
        Coll = np.dot(A_Coll, Coll)
        Coll = coll_reac(Coll, F, dt, coll_prod, coll_death)
        m_coll = min(Coll)
        coll_norm.append(max(Coll))
        
        # Negative Value Controle
        if (m_c < 0 or m_f < 0):
            print("Erreur : Valeur négative")
            
        Index.append(np.argmax(F)/N)
        
        plot_F_C(F, C, Coll, i, t, max_val)  # Plot current iteration
    
    plt.show()
    
# %% Plot
plt.figure()
plt.plot(t[:len(fib_norm)], fib_norm, label=" Fibro. ")
plt.plot(t[:len(cytok_norm)], cytok_norm, label="Cytok.")
plt.plot(t[:len(coll_norm)], coll_norm, label="Coll.")
plt.xlabel("Temps (j)")
plt.ylabel("Max/jour")
plt.title("Val. max Fibro. & Cyto. & Coll.")
plt.legend()
plt.grid()
plt.show()

'''
plt.plot(t,Index)
plt.title("Emplacement du max de fibro.")
plt.xlabel("Temps (j)")
plt.ylabel("x")
plt.show()
'''

    

    
