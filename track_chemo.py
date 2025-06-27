import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Nx, Ny, Nz = 3,3,3
T = 20

liste1 = [np.random.rand(Nx, Ny, Nz) for _ in range(T)]
liste2 = [np.random.rand(Nx, Ny, Nz) for _ in range(T)]

traj_x = []
traj_y = []
traj_z = []

for i in range(T):
    max_mat = np.maximum(liste1[i], liste2[i])
    pos_max = np.unravel_index(np.argmax(max_mat), max_mat.shape)
    traj_x.append(pos_max[0])
    traj_y.append(pos_max[1])
    traj_z.append(pos_max[2])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Création d'un colormap basé sur le temps
colors = cm.viridis(np.linspace(0, 1, T))

for i in range(T-1):
    ax.plot(traj_x[i:i+2], traj_y[i:i+2], traj_z[i:i+2], color=colors[i])

# Scatter pour les points avec couleur selon le temps
p = ax.scatter(traj_x, traj_y, traj_z, c=np.arange(T), cmap='viridis')

fig.colorbar(p, ax=ax, label='Temps')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Trajectoire 3D du point max colorée selon le temps')

plt.show()



import numpy as np

def track_chemotaxis_trajectory(matrix_list1, matrix_list2):
    """
    Suivre la trajectoire 3D du point de valeur maximale issue du max élément par élément
    entre deux listes de matrices 3D, au fil du temps.

    Args:
        matrix_list1 (list of np.ndarray): liste de matrices 3D (Nx x Ny x Nz)
        matrix_list2 (list of np.ndarray): même taille que matrix_list1

    Returns:
        traj_x, traj_y, traj_z (list of int): coordonnées du point max à chaque instant
    """

    assert len(matrix_list1) == len(matrix_list2), "Les listes doivent avoir la même longueur"

    traj_x = []
    traj_y = []
    traj_z = []

    for i in range(len(matrix_list1)):
        max_mat = np.maximum(matrix_list1[i], matrix_list2[i])
        pos_max = np.unravel_index(np.argmax(max_mat), max_mat.shape)
        traj_x.append(pos_max[0])
        traj_y.append(pos_max[1])
        traj_z.append(pos_max[2])

    return traj_x, traj_y, traj_z
    