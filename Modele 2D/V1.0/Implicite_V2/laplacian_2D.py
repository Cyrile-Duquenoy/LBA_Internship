import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal

def Laplacian_2d(N,lbda):
    size = N * N
    diag_main = np.ones(size) * (1 + 4 * lbda)
    diag_x = np.ones(size - 1) * -lbda
    diag_y = np.ones(size - N) * -lbda
    
    for i in range(1, N):
        diag_x[i * N - 1] = 0  # Gestion des bords

    A = sp.diags([diag_main, diag_x, diag_x, diag_y, diag_y], [0, -1, 1, -N, N], format='csr')
    
    # Neumann Conditions
    for i in range(N):  
        for j in range(N):
            index = i * N + j  # Index global dans le vecteur aplati
            
            if i == 0:  # Bord bas
                A[index, index + N] += lbda  # Copie du point adjacent
            if i == N - 1:  # Bord haut
                A[index, index - N] += lbda
            if j == 0:  # Bord gauche
                A[index, index + 1] += lbda
            if j == N - 1:  # Bord droit
                A[index, index - 1] = lbda
     
    return A

