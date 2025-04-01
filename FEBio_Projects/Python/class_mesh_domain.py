class mesh_domain:
    def __init__(self, name_: str, mat_: str):
        self.name_ = name_
        self.mat_ = mat_
        
    def get_name(self):
        return self.name_
    
    def get_mat(self):
        return self.mat_
    
    def set_name(self, name: str):
        self.name_ = name
        
    def set_mat(self, mat: str):
        self.mat_ = mat
    
if __name__ == "__main__":
    M1 = mesh_domain("Part1", "Material1")
    print(M1.get_name())
