import numpy as np
import matplotlib.pyplot as plt
from mesh import Mesh
import os
from pathlib import Path

from febio_python.feb import Feb
from febio_python.feb import Feb40
from febio_python.core import (
    Nodes,
    Elements,
    ShellDomain,
    NodalLoad,
    LoadCurve,
    NodeSet,
    Material,
    SolidDomain
)
from febio_python.feb import run

def data_dict(filename):
    data_dict = {}

    with open(filename, "r") as f:
        for line in f:
            if line.strip():  # ignore les lignes vides
                key, value = line.strip().split()
                data_dict[int(key)] = float(value)
    
    return data_dict

def compute_young_list(coll_dict:dict):
    coll_max = max(coll_dict.values())
    young_min = coll_max / 10
    Y = [float(young_min * n) for n in range(11)]
    return Y
    
def compute_elements_from_young(feb, young_list:list, Coll_data:dict):
    E_list = []
    for ids in range(len(young_list) - 1):
        E = []
        for elmt in Coll_data:
            if young_list[ids] <= Coll_data[elmt] + 30 <= young_list[ids + 1]:
                E.append(elmt)
        E_list.append(E)
    return E_list

                
#def compute_elements_from_young(feb, young_list:list, Coll_data:dict):
    '''
    Les eléments E_i ont un module de Young Y_i qui dépendent d eleur Part_i
    Si le module de Young de E_i + celui de Coll_i dépasse un certain seuil,
    Alors E_i pascule en Part_{i+1}
    '''


if __name__ == "__main__":
    mesh = Mesh('cube.msh')
    elements = mesh.get_elements()
    #print(elements)
    
    Y_init = 30
    Coll_init = 0
    
    feb_path = Path(__file__).parent / 'virgin_iso_elastic.feb'
    feb = Feb(filepath=feb_path)
    print(feb)
    
    
    Coll_interp = data_dict("coll_value_interpolate.txt")
    
    max_Coll = max(Coll_interp.values())
    print(max_Coll)
    
    young_list = compute_young_list(Coll_interp)
    print(young_list)
    print('\n')
    
    
    # Initialisation du dictionnaire final
    interval_dicts = [{} for _ in range(len(young_list)-1)]
    
    # Construction du nouveau format : {moyenne: [éléments]}
    for i in range(len(young_list)-1):
        a = young_list[i]
        b = young_list[i+1]
        mean = (a + b) / 2 
        interval_dicts[i][mean] = []
    
    # Remplissage
    for key, value in Coll_interp.items():
        for i in range(len(young_list)-1):
            a = young_list[i]
            b = young_list[i+1]
            if a <= value < b:
                mean = (a + b) / 2 
                interval_dicts[i][mean].append(key)
                break
    
    # Résultat
    print(interval_dicts)

    
        