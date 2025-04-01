from __future__ import annotations  # Pour éviter les problèmes de référence circulaire
from class_element import element

class part:
    def __init__(self, number_: int):
        self.number_ = number_
        self.elements_ = {}  # Dictionnaire {id: element}

    def get_number(self)->int:
        return self.number_

    def set_number(self, i: int):
        self.number_ = i

    def get_nb_elements(self)->int:
        return len(self.elements_)

    def add_element(self, element_id: int, new_element: element):
        if element_id in self.elements_:
            raise ValueError(f"L'élément avec l'ID {element_id} existe déjà.")
        self.elements_[element_id] = new_element

    def remove_element(self, element_id: int):
        if element_id in self.elements_:
            del self.elements_[element_id]
        else:
            raise KeyError(f"L'élément avec l'ID {element_id} n'existe pas.")

    def get_element(self, element_id: int)->int:
        return self.elements_.get(element_id, None)
        
if __name__ == "__main__":
    e1 = element(10, [1,2,3,4,5,6,7,8], 8, 1)    
    x = e1.get_nodes()
    print(x)
    
    P1 = part(1)
    P1.add_element(1, e1)
    print(P1.get_element(1).get_nodes())
s