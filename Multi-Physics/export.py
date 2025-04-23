import numpy as np
import vtk
from febio_python.feb import Feb
from febio_python.feb import run
from febio_python.core import NodalLoad, LoadCurve, NodeSet
from febio_python.feb import Feb40
from febio_python.core import Nodes, Elements
from febio_python.core import ShellDomain

class Export:
    """
    Classe pour charger et traiter les données de tenseurs à partir de fichiers VTK et FEB.
    """
    def __init__(self, vtk_filename, feb_filename):
        """
        Initialise l'exportateur de données.
        
        Args:
            vtk_filename: Nom du fichier VTK contenant les tenseurs
            feb_filename: Nom du fichier FEB associé
        """
        # === Charger le fichier .vtk ===
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_filename)
        reader.Update()
        
        # === Accès aux données du maillage ===
        self.data = reader.GetOutput()
        self.nb_points = self.data.GetNumberOfPoints()
        self.vtk_filename = vtk_filename
        self.feb_filename = feb_filename
        
        # Charger le fichier FEB si nécessaire
        try:
            from febio_python.feb import Feb
            self.feb = Feb(feb_filename)
        except ImportError:
            print("Warning: Module febio_python non disponible. Certaines fonctionnalités seront limitées.")
            self.feb = None
        except Exception as e:
            print(f"Erreur lors du chargement du fichier FEB {feb_filename}: {e}")
            self.feb = None
        
    def get_nodes_coordinates(self):
        """
        Récupère les coordonnées des nœuds depuis le fichier FEB.
        
        Returns:
            Liste des coordonnées des nœuds [[x1,y1,z1], [x2,y2,z2], ...]
        """
        if self.feb is None:
            # Si le fichier FEB n'est pas disponible, essayer d'extraire les coordonnées du VTK
            coords = []
            for i in range(self.nb_points):
                coords.append(self.data.GetPoint(i))
            return coords
        
        try:
            E = self.feb.get_nodes()
            return E[0].coordinates
        except Exception as e:
            print(f"Erreur lors de la récupération des coordonnées: {e}")
            # Retour de secours: extraire les coordonnées du VTK
            coords = []
            for i in range(self.nb_points):
                coords.append(self.data.GetPoint(i))
            return coords
        
    def get_nb_point(self):
        """
        Récupère le nombre de points.
        
        Returns:
            Nombre de points dans le maillage
        """
        return self.nb_points
    
    def set_nb_point(self):
        """
        Met à jour le nombre de points depuis les données VTK.
        """
        self.nb_points = self.data.GetNumberOfPoints()
        
    def get_tensors(self):
        """
        Récupère les tenseurs des points.
        
        Returns:
            Liste des tenseurs sous forme de matrices 3x3
        """
        tensors = []
        tensor_array = self.data.GetPointData().GetTensors()
        
        if tensor_array is None:
            # Si aucun tenseur n'est trouvé, essayer avec d'autres noms
            for array_name in self.data.GetPointData().GetArrayNames():
                if self.data.GetPointData().GetArray(array_name).GetNumberOfComponents() == 9:
                    tensor_array = self.data.GetPointData().GetArray(array_name)
                    print(f"Utilisation de l'array '{array_name}' comme tenseur")
                    break
        
        if tensor_array is None:
            raise ValueError("Aucun tenseur trouvé dans les données VTK.")
            
        for i in range(self.nb_points):
            flat_tensor = [tensor_array.GetComponent(i, j) for j in range(9)]
            matrix = np.array(flat_tensor).reshape((3, 3))
            tensors.append(matrix)
            
        return tensors
    
    def compute_tensor_norms(self, tensors):
        """
        Calcule les normes des tenseurs.
        
        Args:
            tensors: Liste des tenseurs (matrices 3x3)
            
        Returns:
            Liste des normes des tenseurs
        """
        norms = []
        for tensor in tensors:
            norm = np.sqrt(np.sum(tensor**2))  # norme de Frobenius
            norms.append(norm)
        return norms
    
    def get_principal_stresses(self, tensors):
        """
        Calcule les contraintes principales pour chaque tenseur.
        
        Args:
            tensors: Liste des tenseurs de contrainte
            
        Returns:
            Liste des contraintes principales pour chaque tenseur
        """
        principal_stresses = []
        for tensor in tensors:
            # Calcul des valeurs propres (contraintes principales)
            eigenvalues, _ = np.linalg.eig(tensor)
            # Trier par ordre décroissant
            eigenvalues = np.sort(eigenvalues)[::-1]
            principal_stresses.append(eigenvalues)
        return principal_stresses
    
    def get_von_mises_stresses(self, tensors):
        """
        Calcule les contraintes de von Mises pour chaque tenseur.
        
        Args:
            tensors: Liste des tenseurs de contrainte
            
        Returns:
            Liste des contraintes de von Mises
        """
        von_mises = []
        for tensor in tensors:
            # Extraire les composantes du tenseur
            s11, s22, s33 = tensor[0,0], tensor[1,1], tensor[2,2]
            s12, s13, s23 = tensor[0,1], tensor[0,2], tensor[1,2]
            
            # Calcul de la contrainte de von Mises
            vm = np.sqrt(0.5 * ((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 + 
                               6 * (s12**2 + s13**2 + s23**2)))
            von_mises.append(vm)
        return von_mises


if __name__ == "__main__":
    # Exemple d'utilisation
    try:
        exp = Export('true_test.vtk', 'true_test2.feb')
        
        # Récupérer les tenseurs
        tensors = exp.get_tensors()
        print(f"Nombre de tenseurs: {len(tensors)}")
        
        # Calculer les normes
        norms = exp.compute_tensor_norms(tensors)
        print(f"Normes (min, max): ({min(norms):.3e}, {max(norms):.3e})")
        
        # Récupérer les coordonnées des nœuds
        coords = exp.get_nodes_coordinates()
        print(f"Nombre de nœuds: {len(coords)}")
        
    except Exception as e:
        print(f"Erreur: {e}")

