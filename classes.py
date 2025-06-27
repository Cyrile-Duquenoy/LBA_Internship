class Material:
    def __init__(self, ids:int = None, name:str = None, types:str = 'isotropic elastic', density=1, v=0.3, E = None):
        self.ids = ids
        self.name = name
        self.types = types
        self.density = density
        self.v = v
        self.E = E
        self.attributes = {'E': self.E, 'v': self.v, 'density': self.density}
        
    def __str__(self):
        return f'Material(ids={self.ids}, name={self.name}, types={self.types}, attributes={self.attributes})'
    
class Element:
    def __init__(self,elmt_data, types='hex8',  name: str = None, part: str = None, mat:str = None):
        self.type = types
        self.name = name
        self.elmt_data = elmt_data
        self.part = part
        self.mat = mat
        
    def __str__(self):
        return f'Element(types={self.type}, name={self.name}, elmt_data={self.elmt_data}, part={self.part}, mat={self.mat})'
    
    
class Node:
    def __init__(self,node_data, name='cube'):
        self.node_data = node_data
        self.name = name
        
    def __str__(self):
        return f'Node(name={self.name}, node_data={self.node_data})'
    
class SolidDomain:
    def __init__(self, name, mat):
        self.name = name
        self.mat = mat
        
    def __str__(self):
        return f'SolidDomain(name={self.name}, mat={self.mat})'

