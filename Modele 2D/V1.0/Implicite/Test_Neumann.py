import numpy as np
import matplotlib.pyplot as plt

N=5
x=np.linspace(0,1,N)
y=x

def Laplacian_2d_V2(N,lbda):
    diag_main = np.ones(N*N) * (1 + 4 * lbda)
    diag_y = np.ones(N*N - N) * -lbda
    diag_x = np.ones(N*N - 1)* -lbda
        
    A = np.diag(diag_main) + np.diag(diag_y, -N) + np.diag(diag_y, N) + np.diag(diag_x, -1) + np.diag(diag_x, 1)
     
    return A


M = Laplacian_2d_V2(N, 1)
print(M)

def f(x,y):
    return (1-x**2)*(1-y**2)

def u_init(N,init):
    U=np.eye(N)
    for i in range(N):
        for j in range(N):
            U[i,j]=f(x[i],y[j])
    return U

U0 = u_init(N,0.06)
print(U0)
U0 = U0.flatten()
print(U0)

U = M @ U0
print(U)

U = U.reshape(N,N)
print(U)

U_new = U
U_new[0, :] = U_new[1, :]
U_new[-1, :] = U_new[-2, :]
U_new[:, 0] = U_new[:, 1]
U_new[:, -1] = U_new[:, -2]

# Coins : moyenne des voisins adjacents
U_new[0, 0] = 0.5 * (U_new[0, 1] + U_new[1, 0])      # Coin bas-gauche
U_new[0, -1] = 0.5 * (U_new[0, -2] + U_new[1, -1])   # Coin haut-gauche
U_new[-1, 0] = 0.5 * (U_new[-1, 1] + U_new[-2, 0])   # Coin bas-droit
U_new[-1, -1] = 0.5 * (U_new[-1, -2] + U_new[-2, -1]) # Coin haut-droit
print(U_new)

