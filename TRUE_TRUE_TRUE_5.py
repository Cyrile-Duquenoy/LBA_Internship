import numpy as np
from mesh import Mesh
from pathlib import Path
from BioToFeb import (data_dict,
                      compute_young_list,
                      compute_elements_from_young,
                      )
from classes import(Material,
                    Element,
                    SolidDomain)

from XML import (insert_element_block_in,
                 insert_material_block_in,
                 insert_solid_domain_in)


#%%

def young_liste(coll:dict):
    c_max = max(coll.values())
    c_min = min(coll.values())
    C = np.linspace(c_min, c_max, 10)
    Y = C 
    return Y

'''
def young_liste(coll:dict, young: list):
    c_max = max(coll.values())
    c_min = min(coll.values())
    C = np.linspace(c_min + min(young), c_max + max(young), 10)
    Y = C 
    return Y
'''

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
    return young

def save_young_list_to(young: list, pathname: Path):
    with open(pathname, "w") as f:
        for y in young:
            f.write(f'{y}\n')
            
def compute_youg_list(old_young_path: Path, coll_interpolate_path: Path):
    old_young = read_young_list_from(old_young_path)
    coll = data_dict(coll_interpolate_path)
    
    for key, value in coll.items():
        coll[key] *= 1000
        
    new_young = young_liste(coll, old_young)
    
    return new_young, coll
            


if __name__ == '__main__':
    mesh = Mesh('cube.msh')
    #coll_path = Path(__file__).parent / 'coll_value_interpolate.txt'
    coll_path = 'Z:/Model/coll_value_interpolate_1.txt'
    
    nodes = mesh.get_nodes()
    elements = mesh.get_elements()
    
    '''Récupération des collagènes interpolées sous la forme d'un dictionnaire'''
    coll = data_dict(coll_path)
    for key, value in coll.items():
        if value < 0.050 :
            coll[key] = value * 1000 + 30
        else :
            coll[key] = 80
    print(coll)
    print('\n')
    
    '''Liste de module de Young'''
    young_list = young_list(coll)
    print(young_list)
    print('\n')
    
    '''Initialisation d'une liste de dictionnaire'''
    dico_list = [{} for _ in range(len(young_list)-1)]
    print(dico_list)
    print('\n')
    
    # Construction du nouveau format : {moyenne: [éléments]}
    for i in range(len(young_list)-1):
        a = young_list[i]
        b = young_list[i+1]
        mean = (a + b) / 2
        dico_list[i][mean] = []
        
    print(dico_list)
    
    # Remplissage
    for key, value in coll.items():
        for i in range(len(young_list)-1):
            a = young_list[i]
            b = young_list[i+1]
            if a <= value < b:
                mean = (a + b) / 2
                dico_list[i][mean].append(key)
                if i == 4:
                    print(mean, value)
                break
    print('Interval : \n',dico_list,'\n')
    max_value = max(coll.values())
    if max_value >= young_list[-1]:
        mean = (young_list[-1] + (young_list[-1] + 1)) / 2  # ou tout autre borne supérieure
        dico_list.append({mean: []})
        for key, value in coll.items():
            if value >= young_list[-1]:
                dico_list[-1][mean].append(key)
    print('Interval : \n',dico_list,'\n')
    
    ''''Vérification du nombre d'éléments'''
    for i in range(len(dico_list)):
        for values in dico_list[i].values():
            print('\n')
            print(values)
            print(len(values))
            print('\n')
            
#%%
    feb_path = Path(__file__).parent / 'Z:/Model/iter1.feb'
    
    # Créer le fichier Model.feb si inexistant ou vide
    if not feb_path.exists() or feb_path.stat().st_size == 0:
        root = ET.Element("febio_spec", {"version": "4.0"})
        tree = ET.ElementTree(root)
        tree.write(feb_path, encoding="utf-8", xml_declaration=True)
        print(f"✅ Fichier initialisé : {feb_path}")
    
    
    elmt_data = elmt_data(dico_list)
    
    print('Interval : \n',dico_list,'\n')
    for i in range(len(dico_list)):
        if len(list(dico_list[i].values())[0]) > 0:
            mat = Material(ids=f'{i+1}', name=f'Material{i+1}', E=list(dico_list[i].keys())[0])
            insert_material_block_in(feb_path, mat)
            
            
    
            elmt = Element(elmt_data=elmt_data[i], name=f'Part{i+1}', part=f'Part{i+1}', mat=f'{i+1}')
            insert_element_block_in(feb_path, elmt)
            
            solid_domain = SolidDomain(name=f'Part{i+1}', mat=f'Material{i+1}')
            insert_solid_domain_in(feb_path, solid_domain)
        
    print('\n')
    print(dico_list)
    
        


        

        


            

    

    

            
    

