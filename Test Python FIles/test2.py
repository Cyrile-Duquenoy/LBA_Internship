import numpy as np
import matplotlib.pyplot as pl
import test
from pathlib import Path

def get_new_young_element_set(old_young_element_set, new_young_list):
    if len(old_young_element_set) != len(new_young_list):
        raise ValueError('Les deux listes doivent avoir le même nombre de module de young !')
    new_young_element_set = {}
    for i in range(len(old_young_element_set) - 1):
        a = old_young_element_set[i]
        b = old_young_element_set[i+1]
        mean = (a + b) / 2
    return mean



def get_new_from_init(young_init_path: Path, coll_path: Path):
    old_young_list = read_young_list_from(young_init_path)
    return old_young_list



def get_new(young_path, coll_path, iteration = 'init'):
    if iteration is 'init':
        return get_new_from_init(young_init_path, coll_path)
        
    

if __name__ == '__main__':
    young_init_path = Path('Z:/Automatisation/DoNotTouch_Files/young_init.txt')
    coll_path = coll_path = Path('Z:/Automatisation/results/coll_final_1.txt')
    L = get_new_from_init(young_init_path, coll_path)
            
    

