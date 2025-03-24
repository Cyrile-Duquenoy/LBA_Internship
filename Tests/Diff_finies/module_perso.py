import numpy as np

def norm2(u:list,v:list)->float:
    if(len(u)==len(v)):
        res=0
        for i in range(len(u)):
            res+=u[i]*v[i]
        return np.sqrt(res)


