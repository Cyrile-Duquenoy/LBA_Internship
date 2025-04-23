import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from tqdm import tqdm  # Pour afficher une barre de progression
#from plot_voxels import *

class Diffusion3DSolver:
    """
    Classe pour résoudre un système d'équations de réaction-diffusion en 3D
    avec trois variables: cytokines (C), fibroblastes (F) et collagène (Coll).
    """
    
    def __init__(self, N, L=1.0, T=20, dt=0.001):
        """
        Initialisation du solveur avec les paramètres de simulation.
        
        Parameters:
        -----------
        N : int
            Nombre de points de grille dans chaque dimension
        L : float
            Longueur du domaine cubique
        T : float
            Temps total de simulation (en jours)
        dt : float
            Pas de temps
        """
        # Paramètres de grille
        self.N = N
        self.L = L
        self.dx = L / (N - 1)
        
        # Paramètres temporels
        self.T = T
        self.dt = dt
        self.t = np.arange(0, T + dt, dt)
        self.n_steps = len(self.t)
        
        # Grille spatiale
        self.x = np.linspace(0, L, N)
        self.y = np.linspace(0, L, N)
        self.z = np.linspace(0, L, N)
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z)
        
        # Coefficients de diffusion
        self.D_C = 2.58e-2
        self.D_F = 1.46e-7
        self.D_Coll = 4.59e-8
        self.rho = 1e-4
        
        # Calcul des lambdas (stabilité)
        self.lbda_C = self._calculate_lambda(self.D_C, self.dt, self.dx)
        self.lbda_F = self._calculate_lambda(self.D_F, self.dt, self.dx)
        self.lbda_Coll = self._calculate_lambda(self.D_Coll, self.dt, self.dx)
        
        # Paramètres de réaction
        self.fib_val = 1e-2
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
        
        # Variables d'état
        self.C = None  # Cytokines
        self.F = None  # Fibroblastes
        self.Coll = None  # Collagène
        
        # Stockage des résultats
        self.C_norm = []
        self.F_norm = []
        self.Coll_norm = []
        self.C_results = []
        self.F_results = []
        self.Coll_results = []
    
    def _calculate_lambda(self, D, dt, dx):
        """Calcule le paramètre lambda pour la stabilité numérique"""
        return D * self.dt / (self.dx**2)
    
    def is_cfl_satisfied(self):
        """Vérifie si la condition CFL est satisfaite pour la stabilité en 3D"""
        return (self.lbda_C <= 1/6 and 
                self.lbda_F <= 1/6 and 
                self.lbda_Coll <= 1/6)
    
    def report_cfl_status(self):
        """Affiche le statut de la condition CFL"""
        print(f"lambda_C = {self.lbda_C}, lambda_F = {self.lbda_F}, lambda_Coll = {self.lbda_Coll}")
        print("Ces valeurs doivent être inférieures à 1/6 pour garantir la stabilité en 3D.")
    
    def initialize(self):
        """Initialise les champs de concentration"""
        self.C = self._init_cytokines()
        self.F = self._init_fibroblasts()
        self.Coll = self._init_collagen()
        
        # Stockage des valeurs initiales
        self.C_norm = [np.max(self.C)]
        self.F_norm = [np.max(self.F)]
        self.Coll_norm = [np.max(self.Coll)]
        
        self.C_results = [self.C.copy()]
        self.F_results = [self.F.copy()]
        self.Coll_results = [self.Coll.copy()]
    
    def _init_cytokines(self):
        """Initialise le champ de cytokines avec une distribution sphérique"""
        U = np.zeros((self.N, self.N, self.N))
        x0, y0, z0 = self.N // 2, self.N // 2, self.N // 2
        rad = self.N // 2
        
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    if (i - x0)**2 + (j - y0)**2 + (k - z0)**2 < rad**2:
                        U[i, j, k] = self._cyto_func(self.x[i], self.y[j], self.z[k])
        return U
    
    def _init_fibroblasts(self):
        """Initialise le champ de fibroblastes"""
        U = np.zeros((self.N, self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    U[i, j, k] = self._fib_func(self.x[i], self.y[j], self.z[k])
        return U
    
    def _init_collagen(self):
        """Initialise le champ de collagène avec des valeurs nulles"""
        return np.zeros((self.N, self.N, self.N))
    
    def _cyto_func(self, x, y, z):
        """Fonction pour la concentration initiale de cytokines"""
        return (x * (1 - x) + y * (1 - y) + z * (1 - z)) * 1.13e-1
    
    def _fib_func(self, x, y, z):
        """Fonction pour la concentration initiale de fibroblastes"""
        return self.fib_val * x * (1 - y) * (0.5 + 0.5 * z)
    
    def laplacian_3d(self, U):
        """Calcule le laplacien 3D d'un champ avec conditions aux limites de Neumann"""
        L = np.zeros_like(U)
        dx2 = self.dx**2
        
        # Intérieur du domaine
        L[1:-1, 1:-1, 1:-1] = (
            U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
            U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
            U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
            6 * U[1:-1, 1:-1, 1:-1]
        ) / dx2
        
        # Conditions de Neumann
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
        """Applique un pas de diffusion explicite en 3D avec conditions aux limites de Neumann"""
        V = U.copy()
        V[1:-1, 1:-1, 1:-1] = U[1:-1, 1:-1, 1:-1] + lbda * (
            U[:-2, 1:-1, 1:-1] + U[2:, 1:-1, 1:-1] +
            U[1:-1, :-2, 1:-1] + U[1:-1, 2:, 1:-1] +
            U[1:-1, 1:-1, :-2] + U[1:-1, 1:-1, 2:] - 
            6 * U[1:-1, 1:-1, 1:-1]
        )
        
        # Conditions de Neumann
        V[0, :, :] = V[1, :, :]
        V[-1, :, :] = V[-2, :, :]
        V[:, 0, :] = V[:, 1, :]
        V[:, -1, :] = V[:, -2, :]
        V[:, :, 0] = V[:, :, 1]
        V[:, :, -1] = V[:, :, -2]
        
        return V
    
    def cyto_reaction(self, C, F):
        """Terme de réaction pour les cytokines"""
        return self.cyto_prod * F * C - self.cyto_death * C
    
    def fib_reaction(self, F, C, lapC):
        """Terme de réaction pour les fibroblastes"""
        term = self.fib_prod * C * (1 - F / self.sat_F)
        death = self.fib_death * F
        chemotaxis = self.chi * F * lapC
        return term - death - chemotaxis
    
    def coll_reaction(self, Coll, F, lapF):
        """Terme de réaction pour le collagène"""
        term = self.coll_prod * F * (1 - Coll / self.coll_sat)
        death = self.coll_death * Coll
        chemotaxis = self.chi_coll * Coll * lapF
        return term - death - chemotaxis
    
    def check_negative_values(self, step, day):
        """Vérifie et signale les valeurs négatives dans les champs"""
        has_negative = False
        
        if np.min(self.C) < 0:
            has_negative = True
            print(f"Valeurs négatives de C au pas {step}, jour {day}: min = {np.min(self.C)}")
            # Option: corriger les valeurs négatives
            self.C = np.maximum(self.C, 0)
            
        if np.min(self.F) < 0:
            has_negative = True
            print(f"Valeurs négatives de F au pas {step}, jour {day}: min = {np.min(self.F)}")
            self.F = np.maximum(self.F, 0)
            
        if np.min(self.Coll) < 0:
            has_negative = True
            print(f"Valeurs négatives de Coll au pas {step}, jour {day}: min = {np.min(self.Coll)}")
            self.Coll = np.maximum(self.Coll, 0)
            
        return has_negative
    
    def solve(self, store_full_results=False, display_progress=True):
        """
        Résout le système d'équations de réaction-diffusion.
        
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
    
    #def plot_slice(self, step=-1, z_slice=None, save_path=None):
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
        im0 = axes[0].imshow(C[:, :, z_slice], cmap='viridis')
        axes[0].set_title('Cytokines')
        plt.colorbar(im0, ax=axes[0])
        
        im1 = axes[1].imshow(F[:, :, z_slice], cmap='plasma')
        axes[1].set_title('Fibroblastes')
        plt.colorbar(im1, ax=axes[1])
        
        im2 = axes[2].imshow(Coll[:, :, z_slice], cmap='inferno')
        axes[2].set_title('Collagène')
        plt.colorbar(im2, ax=axes[2])
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()

# Exemple d'utilisation
if __name__ == "__main__":
    # Créer le solveur avec les paramètres par défaut
    solver = Diffusion3DSolver(N=11, L=1.0, T=20, dt=0.001)
    
    # Initialiser les champs
    solver.initialize()
    
    # Vérifier et afficher le statut CFL
    if not solver.is_cfl_satisfied():
        solver.report_cfl_status()
        print("La simulation ne peut pas continuer en raison de l'instabilité numérique.")
    else:
        # Résoudre le système
        solver.solve(store_full_results=True)
        
        # Tracer les résultats
        solver.plot_norms()
        
        # Si les résultats complets ont été stockés, on peut tracer une coupe
        if len(solver.C_results) > 1:
            solver.plot_slice()

