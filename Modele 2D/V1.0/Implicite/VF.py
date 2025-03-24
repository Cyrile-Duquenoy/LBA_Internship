import numpy as np
import matplotlib.pyplot as plt

N=20
x=np.linspace(0,1,N)
dx=1/(N-1)


T = 200

dt = 0.00001
t=[]
for i in range(T):
    t.append(i*dt)
    
D_f = 0.1
D_c = 2

lbda_f = D_f*dt/dx**2
lbda_c = D_c*dt/dx**2

def laplacian_1d(N, lbda):
    d = np.ones(N)*(1-2*lbda)
    d1 = np.ones(N-1)*(lbda)
    return np.diag(d) + np.diag(d1, 1) + np.diag(d1, -1)

A = laplacian_1d(N,lbda_f)
print(A)


def fib_init(N, init):
    return np.ones(N)*init

F = fib_init(N, 0.06)

print(F)

F_max = 0.03
d_f = 0.05


for i in range(1, len(t)):
    F = F @ A + F/(F + F_max) - d_f*F
    plt.plot(x,F)
    plt.show()

    

