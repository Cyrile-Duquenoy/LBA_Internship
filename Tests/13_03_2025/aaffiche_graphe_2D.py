import numpy as np
import matplotlib.pyplot as plt


'''
                    #############################
                    ##### AFFICHE GRAPHE 2D #####
                    #############################
'''

N=100
L=1

def f(x,y):
    return x*y

x=np.linspace(0,L,N)

def F(x:np.ndarray)->np.ndarray:
    F=np.diag(np.zeros(len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            F[i,j]=f(i,j)
    return F

F=F(x)

plt.imshow(F, cmap='viridis', origin='lower', extent=[0, L, 0, L])
#plt.colorbar(label="Concentration (cellules/cm²)")
#plt.title(f"Évolution à t = {t[n]:.2f} j")
plt.xlabel("x")
plt.ylabel("y")
plt.pause(0.5)

        

