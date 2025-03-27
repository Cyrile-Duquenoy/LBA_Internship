import numpy as np
import matplotlib.pyplot as plt

N = 10
dx = 1/(N-1)

T = 30
dt = 0.1
t = np.arange(0, T + dt, dt)

D_F = 1.46e-7
D_C = 2.58e-2

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

def cyto(x):
    return x*(1-x)*1.13e-1

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
    return 1e-4
    
def fibro_init(x):
    U = []
    for i in range(len(x)):
        U.append(fibro(x[i]))
    return U
    
F = fibro_init(x)

plt.plot(x,F)
plt.title("Fibro. Initial")
plt.show()

def cyto_reac(C, F, dt, cyto_prod, cyto_death):
    C += dt*C*F*cyto_prod - dt*C*cyto_death
    return C

def laplacian(C, dx):
    U = np.zeros_like(C)
    for i in range(1, len(C)-1):
        U[i] = (U[i-1] - 2*U[i] + U[i+1]) / (dx**2)
        U[0] = U[1]
        U[-1] = U[-2]
    return U

def fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi):
    return F + dt*( fibro_prod*C*(1 - F/sat_f) - chi*laplacian(C,dx) - fibro_death*F )

# %%
    

def plot_F_C(F, C, i, t, max_val):
    plt.clf()  # Effacer la figure précédente
    
    # Déterminer les limites communes
    min_val = 0
    
    if(max(F) <= max_val /2 or max(F) <= max_val /2 ):
        max_val = max_val/2

    plt.subplot(1, 2, 1)
    plt.plot(F, label="Fibroblastes", color='blue')
    plt.plot(C, label="Cytokines", color='red')
    plt.title(f"Fibroblastes & Cytokines (itération {i}/{len(t)-1})")
    plt.ylim(min_val, max_val)  # Appliquer les mêmes limites d'axe
    plt.legend()
    '''
    plt.subplot(1, 2, 2)
    plt.plot(C, label="Cytokines", color='red')
    plt.title(f"Cytokines (itération {i}/{len(t)-1})")
    plt.ylim(min_val, max_val)  # Même échelle que le premier subplot
    plt.legend()
    '''

    plt.pause(0.01)  # Pause pour mise à jour dynamique
    
# %% PARAM.

cyto_prod = 1
cyto_death = 1
fibro_prod = 1
fibro_death = 1
sat_f = 1
chi = 2

cytok_norm = []
fib_norm = []

Index = []

# %% Time Iterations

if is_cfl(lbda_C) and is_cfl(lbda_F): # Verify CFL conditions
    A_C = mat(N, lbda_C)
    A_F = mat(N, lbda_F)
    cytok_norm.append(max(C))
    fib_norm.append(max(F))
    Index.append(np.argmax(F)/N)
    plt.figure(figsize=(10, 5))
    for i in range(1,len(t)):
        C = np.dot(A_C, C)
        C = cyto_reac(C, F, dt, cyto_prod, cyto_death)
        
# %% Repetitive Trauma Simulation

        if i%100 == 0:
            C+=cyto_init(x)
            
# %% 
        
        m_c = min(C)
        
        cytok_norm.append(max(C))
        fib_norm.append(max(F))
    
        F = np.dot(A_F, F)
        F = fibro_reac(F, C, sat_f, dt, fibro_prod, fibro_death, chi)
        m_f = min(F)
        
        # Controle negative value
        if (m_c < 0 or m_f < 0):
            print("Erreur : Valeur négative")
            
        Index.append(np.argmax(F)/N)
    
        plot_F_C(F, C, i, t, 5.5e-2)  # Plot
    
    plt.show()  # Final Plot
    
# %% Plot
    
plt.figure()
plt.plot(t[:len(cytok_norm)], cytok_norm, label="Val. max Cyto. / Jour")
plt.xlabel("Temps (j)")
plt.ylabel("Cytokine/Jour")
#plt.yscale('log')
#plt.xscale('log')
plt.title("Cytokine")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(t[:len(fib_norm)], fib_norm, label="Val. max Fibro. / Jour")
plt.xlabel("Temps (j)")
plt.ylabel("Fibroblaste/Jour")
#plt.yscale('log')
#plt.xscale('log')
plt.title("Fibroblaste")
plt.legend()
plt.grid()
plt.show()

plt.plot(t,Index)
plt.title("Emplacement du max de fibro.")
plt.xlabel("Temps (j)")
plt.ylabel("x")
plt.show()

    

    
