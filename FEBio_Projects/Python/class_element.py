import numpy as np
import matplotlib.pyplot as plt

class element:
    def __init__(self, id_: int, nodes_: list, nb_node_: int, part_: int):
        self.id_ = id_
        self.nodes_ = nodes_
        self.nb_node_ = nb_node_
        self.part_ = part_
        
    def get_id(self)->int:
        return self.id_
    
    def set_id(self, i: int):
        self.id_ = i
        
    def get_nb_node(self):
         return self.nb_node_
        
    def get_nodes(self):
        if(self.get_nb_node() == len(self.nodes_)):
            return self.nodes_
        else : 
            print("Number of node different of the length of the element")
    
    def set_nodes(self, t:list):
        if len(t) == len(self.nodes_):
            for i in range(len(self.nodes_)):
                self.nodes_[i] = t[i]
                
    def get_part(self)->int:
        return self.part_
        
    def set_part(self, i:int):
        self.part_ = i
    
    
    ''' To write a line as an element in the .feb file '''
    def to_print(self) -> str:
        to_print = f'<elem id="{self.get_id()}">{self.get_nodes()}</elem>'
        to_print = to_print.translate(str.maketrans("[", "]", ""))
        return to_print
        
        
if __name__ == "__main__":
    e1 = element(10, [1,2,3,4,5,6,7,8], 8, 1)    
    x = e1.get_nodes()
    print(x)    
    print(e1.to_print())
    


