import numpy as np
import matplotlib.pyplot as plt
from bio_solve import Diffusion3DSolver
from export import Export

class multi_physics_solve(Diffusion3DSolver, Export):
    def __init__(self, N=11, L=1.0, T=20, dt=0.001, vtk_filename=None, feb_filename=None):
        Diffusion3DSolver.__init__(self, N, L, T, dt)
        
        self.vtk_filename = vtk_filename
        self.feb_filename = feb_filename
        
        if vtk_filename and feb_filename:
            self.exp = Export(vtk_filename, feb_filename)
        else:
            self.exp = None
            
        self.tensors = None
        self.tensor_norms = None
        self.D_C_grid = None
        self.D_F_grid = None
        self.D_Coll_grid = None
        
        if vtk_filename and feb_filename:
            self.load_tensors()
            
        self.lbda_C_grid = None
        self.lbda_F_grid = None
        self.lbda_Coll_grid = None
            
    def load_tensors(self):
        """Charge les tenseurs depuis les fichiers VTK et FEB"""
        try:
            if self.exp is None:
                self.exp = Export(self.vtk_filename, self.feb_filename)
                
            self.tensors = self.exp.get_tensors()
            self.tensor_norms = self.exp.compute_tensor_norms(self.tensors)
            
            print(f"Chargement de {len(self.tensor_norms)} tenseurs réussi")
            
            # Interpoler les normes des tenseurs sur la grille de simulation
            self.interpolate_tensor_norms()
            
        except Exception as e:
            print(f"Erreur lors du chargement des tenseurs: {e}")
            # Utiliser une distribution uniforme par défaut
            self.use_uniform_diffusion()
            
    def interpolate_tensor_norms(self):
        """
        Interpole les normes des tenseurs sur la grille de simulation.
        Si les tenseurs ne sont pas disponibles ou l'interpolation échoue,
        utilise une diffusion uniforme.
        """
        try:
            if self.tensor_norms is None or len(self.tensor_norms) == 0:
                raise ValueError("Aucun tenseur disponible pour l'interpolation")
            
            node_coordinates = self.exp.get_nodes_coordinates()
            
            if len(node_coordinates) != len(self.tensor_norms):
                raise ValueError(f"Incohérence entre les coordonnées ({len(node_coordinates)}) et les tenseurs ({len(self.tensor_norms)})")
            
            # Normaliser les normes des tenseurs pour avoir un facteur d'échelle
            max_norm = max(self.tensor_norms)
            min_norm = min(self.tensor_norms)
            normalized_norms = [(norm - min_norm) / (max_norm - min_norm) if max_norm > min_norm else 0.5 
                               for norm in self.tensor_norms]
            
            # Créer des grilles de coefficients de diffusion
            self.D_C_grid = np.zeros((self.N, self.N, self.N))
            self.D_F_grid = np.zeros((self.N, self.N, self.N))
            self.D_Coll_grid = np.zeros((self.N, self.N, self.N))
            
            # Calculer les coordonnées physiques des points de la grille
            dx = self.L / (self.N - 1)
            grid_points = np.array([[[
                [i * dx, j * dx, k * dx]
                for k in range(self.N)]
                for j in range(self.N)]
                for i in range(self.N)])
            
            # Pour chaque point de la grille, trouver le nœud le plus proche
            for i in range(self.N):
                for j in range(self.N):
                    for k in range(self.N):
                        grid_point = grid_points[i, j, k]
                        
                        # Calculer les distances à tous les nœuds
                        distances = [np.linalg.norm(np.array(grid_point) - np.array(node)) 
                                    for node in node_coordinates]
                        
                        # Trouver l'indice du nœud le plus proche
                        nearest_idx = np.argmin(distances)
                        
                        # Appliquer le facteur d'échelle basé sur la norme du tenseur
                        factor = 0.5 + normalized_norms[nearest_idx]  # Entre 0.5 et 1.5
                        
                        # Modifier les coefficients de diffusion localement
                        self.D_C_grid[i, j, k] = self.D_C * factor
                        self.D_F_grid[i, j, k] = self.D_F * factor
                        self.D_Coll_grid[i, j, k] = self.D_Coll * factor
                        
            print("Interpolation des tenseurs sur la grille réussie")
            
        except Exception as e:
            print(f"Erreur lors de l'interpolation des tenseurs: {e}")
            self.use_uniform_diffusion()
    
    def use_uniform_diffusion(self):
        """
        Configure une diffusion uniforme sur tout le domaine en utilisant
        les coefficients de diffusion par défaut.
        """
        print("Utilisation d'une diffusion uniforme par défaut")
        
        # Créer des grilles uniformes avec les coefficients par défaut
        self.D_C_grid = np.ones((self.N, self.N, self.N)) * self.D_C
        self.D_F_grid = np.ones((self.N, self.N, self.N)) * self.D_F
        self.D_Coll_grid = np.ones((self.N, self.N, self.N)) * self.D_Coll
        
        
    def _calculate_lambda_grid(self, D_grid, dt, dx):
        lbda_grid = np.zeros((self.N, self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    lbda_grid[i, j, k] = self._calculate_lambda(D_grid[i, j, k], dt, dx)
        return lbda_grid
    
    def is_clf(self):
        maxC = np.max(self.lbda_C_grid)
        maxF = np.max(self.lbda_F_grid)
        maxColl = np.max(self.lbda_Coll_grid)
        return (maxC <= 1/6 and maxF and 1/6 and maxColl <= 1/6)
            
    def report_cfl(self):
        """Affiche le statut de la condition CFL"""
        print("Plus grandes valeurs de lambda : ")
        print(f"lambda_C = {np.max(self.lbda_C_grid)}, lambda_F = {np.max(self.lbda_F_grid)}, lambda_Coll = {np.max(self.lbda_Coll_grid)}")
        print("Ces valeurs doivent être inférieures à 1/6 pour garantir la stabilité en 3D.")
    
    
    
    # Fonction Solve à FAIRE !!!!!
    '''
    def solve_with_variable_diffusion(self, store_full_results=False, display_progress=True):
        """
        Résout l'équation de diffusion avec des coefficients variables
        basés sur les tenseurs interpolés.
        """
        # S'assurer que les grilles de diffusion sont initialisées
        if self.D_C_grid :
            self.use_uniform_diffusion()
        
        # Résoudre l'équation de diffusion en utilisant les grilles de diffusion variables
        # Cette méthode devrait remplacer ou étendre la méthode solve() de Diffusion3DSolver
        # avec un traitement tenant compte des coefficients variables
        if not self.is_cfl():
            print("Erreur: La condition CFL n'est pas satisfaite.")
            self.report_cfl()
            return False
        
        print("Début de la simulation...")
        
        iterator = tqdm(range(1, self.n_steps)) if display_progress else range(1, self.n_steps)
        
        for i in iterator:
            day = i * self.dt
            
            # Diffusion et réaction pour C (cytokines)
            self.C = self.diffuse_explicit_3d(self.C, self.lbda_C)
            lapC = self.laplacian_3d(self.C)
            self.C += self.dt * self.cyto_reaction(self.C, self.F)
            
            # Diffusion et réaction pour F (fibroblastes)
            self.F = self.diffuse_explicit_3d(self.F, self.lbda_F)
            lapF = self.laplacian_3d(self.F)
            self.F += self.dt * self.fib_reaction(self.F, self.C, lapC)
            
            # Diffusion et réaction pour Coll (collagène)
            self.Coll = self.diffuse_explicit_3d(self.Coll, self.lbda_Coll)
            self.Coll += self.dt * self.coll_reaction(self.Coll, self.F, lapF)
            
            # Vérifier les valeurs négatives et corriger si nécessaire
            self.check_negative_values(i, day)
            
            # Enregistrer les normes
            self.C_norm.append(np.max(self.C))
            self.F_norm.append(np.max(self.F))
            self.Coll_norm.append(np.max(self.Coll))
            
            # Stockage optionnel des champs complets
            if store_full_results:
                self.C_results.append(self.C.copy())
                self.F_results.append(self.F.copy())
                self.Coll_results.append(self.Coll.copy())
        
        print("Simulation terminée avec succès.")
        return True
        
        # Implémentation de la résolution avec diffusion variable...
        # (Cette partie dépend de l'implémentation de Diffusion3DSolver)
        
        print("Résolution avec diffusion variable effectuée")
        '''
        
    
    def visualize_diffusion_coefficients(self, plane='xy', slice_idx=None):
        """
        Visualise les coefficients de diffusion sur une coupe 2D.
        
        Args:
            plane: Plan de coupe ('xy', 'xz' ou 'yz')
            slice_idx: Indice de la coupe (si None, utilise le milieu)
        """
        if self.D_C_grid is None:
            print("Les coefficients de diffusion ne sont pas initialisés")
            return
            
        if slice_idx is None:
            slice_idx = self.N // 2
            
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['D_C', 'D_F', 'D_Coll']
        grids = [self.D_C_grid, self.D_F_grid, self.D_Coll_grid]
        
        for i, (ax, title, grid) in enumerate(zip(axes, titles, grids)):
            if plane == 'xy':
                im = ax.imshow(grid[:, :, slice_idx], origin='lower')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_title(f'{title} (Z={slice_idx})')
            elif plane == 'xz':
                im = ax.imshow(grid[:, slice_idx, :], origin='lower')
                ax.set_xlabel('X')
                ax.set_ylabel('Z')
                ax.set_title(f'{title} (Y={slice_idx})')
            elif plane == 'yz':
                im = ax.imshow(grid[slice_idx, :, :], origin='lower')
                ax.set_xlabel('Y')
                ax.set_ylabel('Z')
                ax.set_title(f'{title} (X={slice_idx})')
                
            plt.colorbar(im, ax=ax)
            
        plt.tight_layout()
        plt.show()
        
    def export_results_with_tensors(self, filename):
        """
        Exporte les résultats de la simulation avec les informations des tenseurs.
        
        Args:
            filename: Nom du fichier de sortie
        """
        # Cette méthode pourrait combiner les résultats de la simulation avec 
        # les données des tenseurs pour une analyse plus complète
        
        # Implémentation de l'export...
        print(f"Résultats exportés vers {filename}")


if __name__ == "__main__":
    simulation = multi_physics_solve(
        N=10,                      # Taille de la grille 3D
        L=1.0,                     # Taille physique du domaine
        T=20,                      # Temps total (jours)
        dt=0.001,                  # Pas de temps
        vtk_filename='true_test.vtk',     # Fichier VTK avec les tenseurs
        feb_filename='true_test2.feb'     # Fichier FEB associé
    )
    
    # Visualiser les coefficients de diffusion
    simulation.visualize_diffusion_coefficients()
    
    # Résoudre avec diffusion variable
    simulation.solve_with_variable_diffusion()
    
    # Exporter les résultats
    simulation.export_results_with_tensors('resultats_simulation.vtk')
    
    print(simulation._calculate_lambda_grid(np.ones((11,11,11)), simulation.dt, simulation.dx))
    
        
    

