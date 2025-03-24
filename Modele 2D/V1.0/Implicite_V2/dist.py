import numpy as np

'''
arr = np.array([10, 25, 8, 42, 17])
indice_max = np.argmax(arr)

print(indice_max)
'''

def ind_max(A,N):
    B = A.reshape(N,N)
    ind_max_2d = np.unravel_index(np.argmax(B), B.shape)
    return ind_max_2d






'''

arr_2d = np.array([[10, 25, 8], [42, 17, 33]])
indice_max_2d = np.unravel_index(np.argmax(arr_2d), arr_2d.shape)

print(indice_max_2d)  # Affiche (1, 0) pour la valeur 42
'''