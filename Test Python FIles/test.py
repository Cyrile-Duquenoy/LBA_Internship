import numpy as np
import matplotlib.pyplot

from pathlib import Path 
from classes import *
from mesh import Mesh
from BioToFeb import data_dict
from cyto_to_text import *
from coll_interpolate import *

def elmt_dta(dico_list, mesh_path):
    mesh = Mesh(mesh_path)
    L = []
    for i in range(len(dico_list)):
        dico = {}
        for sublist in dico_list[i].values():  # chaque valeur est une liste
            for elmt_id in sublist:  # on boucle sur les identifiants
                dico[elmt_id] = mesh.get_elements()[elmt_id]
        L.append(dico)
    return L

def read_young_list_from(pathname: Path):
    young = []
    with open(pathname, "r") as f:
        for line in f:
            y = float(line.strip())
            young.append(y)
    return np.array(young)

def save_young_list_to(young: np.ndarray, pathname: Path):
    with open(pathname, "w") as f:
        for y in young:
            f.write(f'{y}\n')
            
def data_dict(pathname):
    data_dict = {}
    with open(pathname, "r") as f:
        for line in f:
            if line.strip():  # ignore les lignes vides
                key, value = line.strip().split()
                data_dict[int(key)] = float(value)
    return data_dict




def young_liste(coll:dict, young: list):
    c_max = max(coll.values())
    c_min = min(coll.values())
    C = np.linspace(c_min + min(young), c_max + max(young), len(young) + 1)
    Y = C 
    return Y

def compute_young_list(coll_interpolate: dict, old_young_path: Path, coeff: float=None):
    old_young = read_young_list_from(old_young_path)
    N = len(old_young) 
    if coeff is None: coeff = 1000
    coll_min = min(coll_interpolate.values()) * coeff + old_young[0]
    coll_max = max(coll_interpolate.values()) * coeff + old_young[-1]
    return np.linspace(coll_min, coll_max, N)

def coll_for_young(coll_interpolate: dict, coeff: float=None) ->dict:
    if coeff is None: coeff = 1000
    coll_for_young = dict({})
    for key, value in coll_interpolate.items():
        coll_for_young[key] = value*coeff
    return coll_for_young
    
    
    
'''
def new_dict(new_young, coll_interpolate) -> dict:
    dico_list = [{young:[]}for young in new_young]
    for i in range(len(new_young) - 1):
        a = new_young[i]
        b = new_young[i+1]
        mean = (a + b) / 2
        coll2young = coll_for_young(coll_interpolate)
        for key, value in coll2young.items():
            if a <= value + mean <= b:
                dico_list[i][new_young].append(value)
    return dico_list
'''

def element_set_list(coll_interpolate: dict, new_young_list):
    elmt_set_list = []
    for key, value in coll_interpolate.items():
        for i in range(len(new_young_list) - 1):
            set_list = []



if __name__ == '__main__':
    
    mesh_path = Path('Z:/Automatisation/DoNotTouch_Files/cube.msh')
    coll_path = Path('Z:/Automatisation/results/coll_final_1.txt')
    coll_interpolate_path = Path('Z:/Automatisation/results/coll_value_interpolate1.txt')
    old_young_path = Path('Z:/Automatisation/DoNotTouch_Files/young_init.txt')
    young_init_path = old_young_path
    new_young_path = Path('Z:/Automatisation/results/new_young.txt')
    
    coll = load_from_file(coll_path)
    print(coll)
    print('\n')
    
    coll_interpolate = compute_coll_interpolate(mesh_path, coll_path)
    print(coll_interpolate)
    print('\n')
    
    new_young = compute_young_list(coll_interpolate, old_young_path, coeff=0)
    save_young_list_to(new_young, new_young_path)
    print(new_young)
    print('\n')
    
    coll2young = coll_for_young(coll_interpolate)
    print(coll)
    
    elements = Mesh(mesh_path).get_elements()
    print(elements)
    
    nodes = Mesh(mesh_path).get_nodes()
    print(nodes)
    
    print(coll_interpolate)
    
    print(new_young)
    

    
    

    

    
    
    
    
    
    

    
    
    
    

