import numpy as np
import matplotlib.pyplot as plt

D = 0.01 #Coeff. Diffusion

L = 1 # Longueur du domaine
N = 11 # Nombres de points de discretisation en espace

T = 200.  #Temps Max

x = np.linspace(0,L,N) # Discretisation du domaine
dx = L/(N-1)

dt=0.5

lbda=D*dt/(dx**2)

eps=10e-6 #Critère d'arret

def mat(N:int, lbda:float)->list:
    d=np.ones(N)*(1+2*lbda)
    d1=np.ones(N-1)*(-lbda)
    return np.diag(d)+np.diag(d1,-1)+np.diag(d1,1)

def f(x:float)->float:
    return x*(1-x)

def F(x:list)->list:
    FF=[]
    for i in range(len(x)):
        FF.append(f(x[i]))
    return FF

F=F(x)
plt.plot(x,F)
plt.title("Test F")
plt.show()

U=np.ones(N)

def norm_L2(u:list,v:list)->float:
    if(len(u)==len(v)):
        res=0
        for i in range(len(u)):
            res+=u[i]*v[i]
        return np.sqrt(res)
    
A=mat(N,lbda)
t=0
stop=eps
UU=[U]
NN=[norm_L2(U,U)]
tt=[0]
while(t<=T or stop<eps):
    count=1
    U=np.linalg.solve(A,U)
    t=t+dt
    plt.plot(x,U)
    plt.show()
    stop=norm_L2(U,UU[count-1])
    NN.append(norm_L2(U,U))
    tt.append(t)
    t+=dt
    
plt.plot(tt,NN)
plt.title("Norme L2")
plt.show()