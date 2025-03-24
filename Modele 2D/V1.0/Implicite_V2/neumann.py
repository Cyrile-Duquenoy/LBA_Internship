import numpy as np

def neumann(U,N):
    U_new = U
    U_new = U_new.reshape(N,N)
    U_new[0, :] = U_new[1, :]
    U_new[-1, :] = U_new[-2, :]
    U_new[:, 0] = U_new[:, 1]
    U_new[:, -1] = U_new[:, -2]

    # Coins : moyenne des voisins adjacents
    U_new[0, 0] = 0.5 * (U_new[0, 1] + U_new[1, 0])      # Coin bas-gauche
    U_new[0, -1] = 0.5 * (U_new[0, -2] + U_new[1, -1])   # Coin haut-gauche
    U_new[-1, 0] = 0.5 * (U_new[-1, 1] + U_new[-2, 0])   # Coin bas-droit
    U_new[-1, -1] = 0.5 * (U_new[-1, -2] + U_new[-2, -1]) # Coin haut-droit
    
    U_new = U_new.flatten()
    
    return U_new
    

