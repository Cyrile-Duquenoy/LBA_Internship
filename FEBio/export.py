import numpy as np
import vtk
from pathlib import Path
from mesh import Mesh

class Export:
    """
    Classe pour charger et traiter les données de tenseurs à partir d'un fichier VTK.
    """
    def __init__(self, vtk_path):
        self.vtk_path = vtk_path

        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(str(vtk_path))  # conversion en str si besoin
        reader.Update()

        self.data = reader.GetOutput()
        self.nb_points = self.data.GetNumberOfPoints()

    def get_nodes_coordinates(self):
        """
        Retourne les coordonnées des nœuds depuis le fichier VTK.
        """
        coords = [self.data.GetPoint(i) for i in range(self.nb_points)]
        return coords

    def get_nb_point(self):
        return self.nb_points

    def set_nb_point(self):
        self.nb_points = self.data.GetNumberOfPoints()

    def get_tensors(self, name="nodal_strain"):
        tensors = []
        point_data = self.data.GetPointData()
        tensor_array = point_data.GetArray(name)

        if tensor_array is None:
            raise ValueError(f"Aucun tableau nommé '{name}' trouvé dans le fichier VTK.")

        if tensor_array.GetNumberOfComponents() != 9:
            raise ValueError(f"Le tableau '{name}' ne contient pas 9 composantes par point (ce n'est donc pas un tenseur 3x3).")

        for i in range(self.nb_points):
            flat_tensor = [tensor_array.GetComponent(i, j) for j in range(9)]
            matrix = np.array(flat_tensor).reshape((3, 3))
            tensors.append(matrix)

        return tensors

    def compute_tensor_norms(self, tensors):
        norms = [np.sqrt(np.sum(tensor**2)) for tensor in tensors]
        #norms = np.linalg.norm(tensors, axis=(1,2))  # shape = (nb_nodes,)
        return norms

    def get_principal_stresses(self, tensors):
        principal_stresses = []
        for tensor in tensors:
            eigenvalues, _ = np.linalg.eig(tensor)
            principal_stresses.append(np.sort(eigenvalues)[::-1])
        return principal_stresses

    def get_von_mises_stresses(self, tensors):
        von_mises = []
        for tensor in tensors:
            s11, s22, s33 = tensor[0, 0], tensor[1, 1], tensor[2, 2]
            s12, s13, s23 = tensor[0, 1], tensor[0, 2], tensor[1, 2]

            vm = np.sqrt(0.5 * ((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 +
                                6 * (s12**2 + s13**2 + s23**2)))
            von_mises.append(vm)
        return von_mises
    

    

if __name__ == '__main__':
    vtk_path = Path('Z:/Model/init.vtk')
    mesh_path = Path('Z:/Model/cube.msh')
    mesh = Mesh(mesh_path)
    exp = Export(vtk_path)
    N = 11
    norms = exp.compute_tensor_norms(exp.get_tensors(name="nodal_strain"))
    nodes = mesh.get_nodes()
    Cyt = np.zeros((N,N,N))
    x = np.linspace(0,1,N)
    y=x
    z=x
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                count = 0
                for value in nodes.values():
                    count +=1 
                    if value == (x[i],y[j],z[k]):
                        Cyt[i,j,k] = norms[count -1]
                        count=0
                        
    print(Cyt)
    print(norms[0])



