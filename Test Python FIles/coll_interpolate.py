import numpy as np
from mesh import Mesh
from cyto_to_text import *

from pathlib import Path

def interpolate(values_list):
    if len(values_list) != 8:
        raise ValueError("Value Error")
    res = 0
    for value in values_list:
        res += value
    return res / 8

def build_interpolated_element_values(Elements, Nodes, Coll):
    """
    Pour chaque élément, interpole la valeur en faisant la moyenne
    des valeurs aux 8 nœuds (définis par leurs coordonnées).
    
    Arguments :
        Elements : dict {elmt_id : list[node_id, ...]}
        Nodes    : dict {node_id : (x, y, z)}
        Coll     : dict {(x, y, z) : valeur}
        
    Retour :
        dict {elmt_id : val_interp}
    """
    Val = {}

    for elmt_id, node_ids in Elements.items():
        values_list = []
    
        for node_id in node_ids:
            coord = Nodes[node_id]
            rounded_coord = tuple(round(c, 6) for c in coord)
    
            if rounded_coord not in Coll:
                raise KeyError(f"La coordonnée {rounded_coord} est absente de Coll.")
            
            values_list.append(Coll[rounded_coord])

        val_interp = interpolate(values_list)
        Val[elmt_id] = val_interp

    return Val


def save_coll_elmt_interpolate(filename, Val):
    """Sauvegarde les valeurs de cytokines avec coordonnées dans un fichier texte."""
    with open(filename, "w") as f:
        for key, value in Val.items():
            f.write(f"{key}\t{value:.6f}\n")
            
def compute_coll_interpolate(mesh_path: Path, coll_path: Path) ->dict:
    coll = load_from_file(coll_path)
    mesh = Mesh(mesh_path)
    elements = mesh.get_elements()
    nodes = mesh.get_nodes()
    coll_interpolate = build_interpolated_element_values(elements, nodes, coll)
    return coll_interpolate
    
            
            
if __name__ == '__main__':

    mesh = Mesh('Z:/Automatisation/DoNotTouch_Files/cube.msh')
    
    
    i = 1
    
    Elements = mesh.get_elements() 
    #print(Elements)
    print('\n')
    
    coll_path = ''
    
    Coll = load_from_file(f'Z:/Model/coll_final_{i}.txt')
    #print(Coll)
    
    
    Nodes = mesh.get_nodes()
    print(Nodes)
    
    
    Val = build_interpolated_element_values(Elements, Nodes, Coll)
    print(Val)
    save_coll_elmt_interpolate('Z:/Model/coll_value_interpolate_1.txt', Val)

