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
     
    return A