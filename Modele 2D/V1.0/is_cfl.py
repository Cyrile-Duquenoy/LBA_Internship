import numpy as np

def is_cfl(dt,dx,D):
    return dt <= D*dx**2/4
    

