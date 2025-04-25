import numpy as np
import matplotlib.pyplot as plt
from export import Export
from tqdm import tqdm

class Liver_Fibrosis_Solver(Export):
    def __init__(self,
                 N=11,
                 L=1.0,
                 dt=0.001,
                 T=20,
                 C=False,
                 F=None,
                 Coll=None,
                 vtk_filename=None,
                 feb_filename=None):
        Export.__init__(self, vtk_filename, feb_filename)
        
        # Mesh Parameters
        self.N = N 
        self.L = L
        self.dx = L / (N -1)
        
        # Time
        self.T = T
        self.dt = dt
        self.t = np.arange(0, T + dt, dt)
        self.n_steps = len(self.t)
        
        # Space
        self.x = np.linspace(0, L, N)
        self.y = np.linspace(0, L, N)
        self.z = np.linspace(0, L, N)
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z)
        
        # Diff. Coeff.
        self.D_C = 2.58e-2
        self.D_F = 1.46e-7
        self.D_Coll = 4.59e-8
        self.rho = 1e-4
        
        # Stability
        self.lbda_C = self._calculate_lambda(self.D_C, self.dt, self.dx)
        self.lbda_F = self._calculate_lambda(self.D_F, self.dt, self.dx)
        self.lbda_Coll = self._calculate_lambda(self.D_Coll, self.dt, self.dx)
        
        # Parameters
        self.fib_val = 1e-2 # Default concentration of fibroblasts
        self.fib_prod = 1
        self.fib_death = 0.1
        self.cyto_prod = 1
        self.cyto_death = 1
        self.coll_prod = 0.6
        self.coll_death = 0.6
        self.chi = 2
        self.sat_F = 1
        self.coll_sat = 1
        self.chi_coll = 1e-12
        
        # Cytokines, Fibroblasts and Collagen
        self.C = C
        self.F = F 
        self.Coll = Coll 
        
        # To stock norms results
        self.C_norm = []
        self.F_norm = []
        self.Coll_norm = []
        
        # To stock results
        self.C_results = []
        self.F_results = []
        self.Coll_results = []
        
        # To use .vtk anf .feb files
        self.export = Export(vtk_filename, feb_filename)

    def load_tensors_norms(self):
        norms = self.export.compute_tensor_norms(self.export.get_tensors())
        return norms
    
    def _calculate_lambda(self, D, dt, dx):
        return D * self.dt / (self.dx**2)
    
    def is_cfl_satisfied(self):
        return (self.lbda_C <= 1/6 and 
                self.lbda_F <= 1/6 and 
                self.lbda_Coll <= 1/6)
    
    def report_cfl_status(self):
        print(f"lambda_C = {self.lbda_C}, lambda_F = {self.lbda_F}, lambda_Coll = {self.lbda_Coll}")
        print("Ces valeurs doivent être inférieures à 1/6 pour garantir la stabilité en 3D.")
        
    def initialize(self):
        if self.C is None:
            self.C = self._init_cytokines()
        if not self.F:
            self.F = self._init_fibroblasts()
        if not self.Coll:
            self.Coll = self._init_collagen()
            
        # Stockage des valeurs initiales
        self.C_norm = [np.max(self.C)]
        self.F_norm = [np.max(self.F)]
        self.Coll_norm = [np.max(self.Coll)]
          
        self.C_results = [self.C.copy()]
        self.F_results = [self.F.copy()]
        self.Coll_results = [self.Coll.copy()]
        
    def _init_cytokines(self):
        if not self.C:
            U = np.zeros((self.N, self.N, self.N))
            x0, y0, z0 = self.N // 2, self.N // 2, self.N // 2
            rad = self.N // 2
            
            for i in range(self.N):
                for j in range(self.N):
                    for k in range(self.N):
                        if (i - x0)**2 + (j - y0)**2 + (k - z0)**2 < rad**2:
                            U[i, j, k] = self._cyto_func(self.x[i], self.y[j], self.z[k])
        return self.C
    
    def _init_fibroblasts(self):
        U = np.zeros((self.N, self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    U[i, j, k] = self._fib_func(self.x[i], self.y[j], self.z[k])
        return U
    
    def _init_collagen(self):
        return np.zeros((self.N, self.N, self.N))
    
    def _cyto_func(self, x, y, z):
        return (x * (1 - x) + y * (1 - y) + z * (1 - z)) * 1.13e-1
    
    def _fib_func(self, x, y, z):
        return self.fib_val
    
    def plot(self, time_index=0):
        if not (self.C_results and self.F_results and self.Coll_results):
            print("Erreur : les résultats doivent être initialisés avec la méthode initialize() avant d'utiliser plot().")
            return
    
        if time_index >= len(self.C_results):
            print(f"Erreur : time_index ({time_index}) hors limites. Maximum autorisé : {len(self.C_results)-1}.")
            return
    
        # S'assurer que les données sont des arrays NumPy
        C = np.array(self.C_results[time_index])
        F = np.array(self.F_results[time_index])
        Coll = np.array(self.Coll_results[time_index])
        center = self.N // 2  # Coupe centrale
    
        fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    
        im0 = axs[0].imshow(C[:, :, center], origin='lower', cmap='viridis')
        axs[0].set_title(f'Cytokines (z={center})')
        fig.colorbar(im0, ax=axs[0])
    
        im1 = axs[1].imshow(F[:, :, center], origin='lower', cmap='plasma')
        axs[1].set_title(f'Fibroblastes (z={center})')
        fig.colorbar(im1, ax=axs[1])
    
        im2 = axs[2].imshow(Coll[:, :, center], origin='lower', cmap='inferno')
        axs[2].set_title(f'Collagène (z={center})')
        fig.colorbar(im2, ax=axs[2])
    
        plt.suptitle(f"Visualisation des champs à t = {self.t[time_index]:.3f}s (index={time_index})", fontsize=14)
        plt.tight_layout()
        plt.show()
        
    def laplacian_3d(self, U):
        L = np.zeros_like(U)
        dx2 = self.dx**2
        
        L[1:-1, 1:-1, 1:-1] = (
            U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
            U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
            U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
            6 * U[1:-1, 1:-1, 1:-1]
        ) / dx2
        
        # Neumann Conditions
        
        # Face x=0
        L[0, :, :] = L[1, :, :]
        # Face x=L
        L[-1, :, :] = L[-2, :, :]
        # Face y=0
        L[:, 0, :] = L[:, 1, :]
        # Face y=L
        L[:, -1, :] = L[:, -2, :]
        # Face z=0
        L[:, :, 0] = L[:, :, 1]
        # Face z=L
        L[:, :, -1] = L[:, :, -2]
        
        return L
    
    def diffuse_explicit_3d(self, U, lbda):
        V = U.copy()
        V[1:-1, 1:-1, 1:-1] = U[1:-1, 1:-1, 1:-1] + lbda * (
            U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
            U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
            U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
            6 * U[1:-1, 1:-1, 1:-1]
        )
        
        # Neumann Conditions
        V[0, :, :] = V[1, :, :]
        V[-1, :, :] = V[-2, :, :]
        V[:, 0, :] = V[:, 1, :]
        V[:, -1, :] = V[:, -2, :]
        V[:, :, 0] = V[:, :, 1]
        V[:, :, -1] = V[:, :, -2]
        
        return V
    
    def cyto_reaction(self, C, F):
        return self.cyto_prod * F * C - self.cyto_death * C
    
    def fib_reaction(self, F, C, lapC):
        term = self.fib_prod * C * (1 - F / self.sat_F)
        death = self.fib_death * F
        chemotaxis = self.chi * F * lapC
        return term - death - chemotaxis
    
    def coll_reaction(self, Coll, F, lapF):
        term = self.coll_prod * F * (1 - Coll / self.coll_sat)
        death = self.coll_death * Coll
        chemotaxis = self.chi_coll * Coll * lapF
        return term - death - chemotaxis
    
    def check_negative_values(self, step, day):
        has_negative = False
        
        if np.min(self.C) < 0:
            has_negative = True
            print(f"Valeurs négatives de C au pas {step}, jour {day}: min = {np.min(self.C)}")
            # Option: corriger les valeurs négatives
            #self.C = np.maximum(self.C, 0)
            
        if np.min(self.F) < 0:
            has_negative = True
            print(f"Valeurs négatives de F au pas {step}, jour {day}: min = {np.min(self.F)}")
            #self.F = np.maximum(self.F, 0)
            
        if np.min(self.Coll) < 0:
            has_negative = True
            print(f"Valeurs négatives de Coll au pas {step}, jour {day}: min = {np.min(self.Coll)}")
            #self.Coll = np.maximum(self.Coll, 0)
            
        return has_negative
        
    def solve(self, store_full_results=False, display_progress=True):
        """
        Résout le système d'équations de réaction-diffusion-advection.
        
        Parameters:
        -----------
        store_full_results : bool
            Si True, stocke les champs complets à chaque pas de temps (utilisation mémoire importante)
        display_progress : bool
            Si True, affiche une barre de progression
        
        Returns:
        --------
        success : bool
            True si la simulation s'est terminée avec succès
        """
        if not self.is_cfl_satisfied():
            print("Erreur: La condition CFL n'est pas satisfaite.")
            self.report_cfl_status()
            return False
        
        print("Début de la simulation...")
        
        # Barre de progression
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
    
    def plot_slice(self, step=-1, z_slice=None, save_path=None):
        """
        Trace une coupe 2D des champs à un pas de temps donné
        
        Parameters:
        -----------
        step : int
            Indice du pas de temps à visualiser (-1 pour le dernier)
        z_slice : int
            Indice de la coupe en z (None pour utiliser N//2)
        save_path : str
            Chemin pour sauvegarder la figure (None pour ne pas sauvegarder)
        """
        if z_slice is None:
            z_slice = self.N // 2
            
            
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        # Assurer que les données sont disponibles
        if len(self.C_results) <= step:
            print(f"Erreur: pas assez de résultats stockés. Utilisez store_full_results=True lors de l'appel à solve().")
            return
            
        # Sélectionner les champs au pas de temps demandé
        C = self.C_results[step]
        F = self.F_results[step]
        Coll = self.Coll_results[step]
        
        # Tracer les coupes
        im0 = axes[0].imshow(C[:, :, z_slice], cmap='inferno')
        axes[0].set_title('Cytokines')
        plt.colorbar(im0, ax=axes[0])
        
        im1 = axes[1].imshow(F[:, :, z_slice], cmap='inferno')
        axes[1].set_title('Fibroblastes')
        plt.colorbar(im1, ax=axes[1])
        
        im2 = axes[2].imshow(Coll[:, :,z_slice], cmap='inferno')
        axes[2].set_title('Collagène')
        plt.colorbar(im2, ax=axes[2])
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
    def plot_norms(self, save_path=None):
        """Trace l'évolution des valeurs maximales au cours du temps"""
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.F_norm, label='Fibroblastes')
        plt.plot(self.t, self.C_norm, label='Cytokines')
        plt.plot(self.t, self.Coll_norm, label='Collagène')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        plt.title("Évolution des concentrations maximales")
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()

#%%
def compute_cytokines(vtk_filename, feb_filename):
    """
    Fonction en dehors de la classe :
    ---------------------------------
    Pour monter la cytokine initiale
    à partir des fichiers .vtk et .feb
    """
    norms = exp.compute_tensor_norms(exp.get_tensors())
    N = 11
    nodes = exp.get_nodes_coordinates()
    Cyt = np.zeros((N, N, N))
    
    for i in range(len(nodes)):
        x, y, z = map(int, nodes[i])
        
        if 0 <= x < N and 0 <= y < N and 0 <= z < N:
            Cyt[x, y, z] = norms[i]
        else:
            print(f"Coordonnées hors bornes ignorées : ({x}, {y}, {z})")
    return Cyt
