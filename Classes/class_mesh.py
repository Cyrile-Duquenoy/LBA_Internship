import numpy as np
from point import *

class mesh:
    def __init__(self, x_min, x_max, y_min, y_max, N, M):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.N = N
        self.M = M
    
    def create_mesh(self):
        return np.diag(np.zeros(self.N))

M = mesh(0,1,0,1,5,5).create_mesh()
print(M)