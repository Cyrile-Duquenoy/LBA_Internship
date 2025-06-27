import numpy as np
import matplotlib.pyplot as plt
from export import Export
from tqdm import tqdm
from cyto_to_text import *
from mesh import Mesh



'''
Simule la biologie à partir
des résultats de la simu méca
qui sont stockés dans le fichier .vtk
'''
class Liver_Fibrosis_Solver(Export):
    def __init__(self,
                 N=11,
                 L=1.0,
                 dt=0.001,
                 T=20,
                 C=False,
                 F=None,
                 Coll=False,
                 vtk_filename=None,
                 feb_filename=None):
        Export.__init__(self, vtk_filename)
        
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
        
        '''
        # Diff. Coeff.
        self.D_C = 2.58e-2
        self.D_F = 1.47e-6
        self.D_Coll = 4.59e-12
        self.rho = 1e-4
        '''
        # Diff. Coeff.
        self.D_C = 1.08e-2
        self.D_F = 1.47e-7
        self.D_Coll = 4.59e-12
        self.rho = 1e-4
        
        # Stability
        self.lbda_C = self._calculate_lambda(self.D_C, self.dt, self.dx)
        self.lbda_F = self._calculate_lambda(self.D_F, self.dt, self.dx)
        self.lbda_Coll = self._calculate_lambda(self.D_Coll, self.dt, self.dx)
        
        '''
        # Parameters
        self.fib_val = 1.2e-2 # Default concentration of fibroblasts
        self.fib_prod = 1
        self.fib_death = 0.1
        self.cyto_prod = 1
        self.cyto_death = 5
        self.coll_prod = 0.6
        self.coll_death = 0.5
        self.chi = 0.001
        self.sat_F = self.fib_prod / 2
        self.coll_sat = self.coll_prod / 2
        self.chi_coll = 1e-2
        '''
        
        # Parameters
        self.fib_val = 4.75e-3 # Default concentration of fibroblasts
        self.coll_val = 3.26e-3
        self.fib_prod = 5e-4
        self.fib_death = 1.66e-2
        self.cyto_prod = 9.64e-2
        self.cyto_death = 2.376e-1
        self.coll_prod = 0.6e-4
        self.coll_death = 0.6e-4
        self.chi = 0.1
        self.sat_F = 3e-2
        self.coll_sat = 10e-2
        self.chi_coll = 1e-2

        
        
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
        self.export = Export(vtk_filename)

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
        if self.F is None:
            self.F = self._init_fibroblasts()
        if self.Coll is None:
            self.Coll = self._init_coll()
            
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
                            self.C = U
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
    
    def _coll_func(self, x, y, z):
        return self.coll_val
    
    def _init_coll(self):
        if not self.Coll:
            U = np.zeros((self.N, self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    for k in range(self.N):
                        U[i, j, k] = self._coll_func(self.x[i], self.y[j], self.z[k])
                        self.Coll = U
        return self.Coll
        
    
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
        prod = self.cyto_prod * C * F / (F + self.sat_F)
        death = self.cyto_death * C
        #return self.cyto_prod * F * C - self.cyto_death * C
        return prod - death
    
    def fib_reaction(self, F, C, lapC):
        '''
        term = self.fib_prod * C * (1 - F / self.sat_F)
        death = self.fib_death * F
        #chemotaxis = self.chi * F * lapC
        chemotaxis = self.chi * (F / (1 + self.sat_F*F)) * lapC
        return term - death - chemotaxie
        '''
        prod = self.fib_prod * C * (1 - F / self.sat_F)
        death = self.fib_death * F
        chemotaxis = self.chi * (F / (1 + self.sat_F * F)) * lapC
        
        return prod - death - chemotaxis

    
    def coll_reaction(self, Coll, F, lapF):
        '''
        #term = self.coll_prod * F * (1 - Coll / self.coll_sat)
        term = self.coll_prod * (F / (self.sat_F + F)) * (Coll / (self.coll_sat + Coll))
        death = self.coll_death * Coll
        #chemotaxis = self.chi_coll * (Coll / (1 -  Coll / self.coll_sat ) ) * lapF
        chemotaxis = self.chi_coll * (Coll / (Coll + self.coll_sat)) * lapF
        return term - death
        '''
        #eps = 1e-6
        #prod = self.coll_prod * F * (Coll / (self.coll_sat + Coll))
        prod = self.coll_prod * F *(1 - Coll / self.coll_sat)
        death = self.coll_death * Coll
        #chemotaxis = self.chi_coll * (Coll / (Coll + self.coll_sat)) * lapF  # évite la division dangereuse
        #$chemotaxis = self.chi_coll * (Coll / (1 + self.coll_sat * Coll)) * lapF
        chemotaxis = 0
        #print(prod - death - chemotaxis)
        return prod - death - chemotaxis
    
        
    
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
        print("\n")
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
            #z_slice = self.N // 2
            z_slice = 0
            
            
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
    

    def plot_slice_custom(self, step=-1, x_slice=None, y_slice=None, z_slice=None, save_path=None):

        # Assurer que les données sont disponibles
        if len(self.C_results) <= step:
            print(f"Erreur: pas assez de résultats stockés. Utilisez store_full_results=True lors de l'appel à solve().")
            return
            
        # Valeurs par défaut pour les coupes
        if x_slice is None:
            x_slice = self.N // 2
        if y_slice is None:
            y_slice = self.N // 2
        if z_slice is None:
            z_slice = self.N // 2
            

        # Sélectionner les champs au pas de temps demandé
        C = self.C_results[step]
        F = self.F_results[step]
        Coll = self.Coll_results[step]

        C_max = np.max([np.max(res) for res in self.C_results])
        F_max = np.max([np.max(res) for res in self.F_results])
        Coll_max = np.max([np.max(res) for res in self.Coll_results])
        
        Coll_min = np.min([np.min(res) for res in self.Coll_results])

        # Créer la figure 3x3
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        
        # Ligne 1: Coupes selon l'axe X (plan YZ)
        im00 = axes[0, 0].imshow(C[x_slice, :, :], cmap='inferno', origin='lower', vmin=0, vmax=C_max)
        axes[0, 0].set_title(f'Cytokines - Coupe X={self.dx*x_slice}')
        plt.colorbar(im00, ax=axes[0, 0])
        
        im01 = axes[0, 1].imshow(F[x_slice, :, :], cmap='inferno', origin='lower', vmin=0, vmax=F_max)
        axes[0, 1].set_title(f'Fibroblastes - Coupe X={self.dx*x_slice}')
        plt.colorbar(im01, ax=axes[0, 1])
        
        im02 = axes[0, 2].imshow(Coll[x_slice, :, :], cmap='inferno', origin='lower', vmax=Coll_max, vmin=Coll_min)
        axes[0, 2].set_title(f'Collagène - Coupe X={self.dx*x_slice}')
        plt.colorbar(im02, ax=axes[0, 2])
        
        # Ligne 2: Coupes selon l'axe Y (plan XZ)
        im10 = axes[1, 0].imshow(C[:, y_slice, :], cmap='inferno', origin='lower', vmin=0, vmax=C_max)
        axes[1, 0].set_title(f'Cytokines - Coupe Y={self.dx*y_slice}')
        plt.colorbar(im10, ax=axes[1, 0])
        
        im11 = axes[1, 1].imshow(F[:, y_slice, :], cmap='inferno', origin='lower', vmin=0, vmax=F_max)
        axes[1, 1].set_title(f'Fibroblastes - Coupe Y={self.dx*y_slice}')
        plt.colorbar(im11, ax=axes[1, 1])
        
        im12 = axes[1, 2].imshow(Coll[:, y_slice, :], cmap='inferno', origin='lower', vmax=Coll_max, vmin=Coll_min)
        axes[1, 2].set_title(f'Collagène - Coupe Y={self.dx*y_slice}')
        plt.colorbar(im12, ax=axes[1, 2])
        
        # Ligne 3: Coupes selon l'axe Z (plan XY)
        im20 = axes[2, 0].imshow(C[:, :, z_slice], cmap='inferno', origin='lower', vmin=0, vmax=C_max)
        axes[2, 0].set_title(f'Cytokines - Coupe Z={self.dx*z_slice}')
        plt.colorbar(im20, ax=axes[2, 0])
        
        im21 = axes[2, 1].imshow(F[:, :, z_slice], cmap='inferno', origin='lower', vmin=0, vmax=F_max)
        axes[2, 1].set_title(f'Fibroblastes - Coupe Z={self.dx*z_slice}')
        plt.colorbar(im21, ax=axes[2, 1])
        
        im22 = axes[2, 2].imshow(Coll[:, :, z_slice], cmap='inferno', origin='lower', vmax=Coll_max, vmin=Coll_min)
        axes[2, 2].set_title(f'Collagène - Coupe Z={self.dx*z_slice}')
        plt.colorbar(im22, ax=axes[2, 2])
        
        # Titre général
        time_value = self.t[step] if step >= 0 else self.t[step]
        plt.suptitle(f"Coupes 2D à t = {time_value:.3f}j", fontsize=16)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()

        
    def plot_norms(self, save_path=None, ids = None):
        """Trace l'évolution des valeurs maximales au cours du temps"""
        plt.rcParams.update({'font.size': 18})
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.F_norm, label='Fibroblastes')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        if ids is None:
            plt.title(f"Évolution des concentrations maximales")
        else :
            plt.title(f"Évolution des concentrations maximales. Iterations {ids}")
        plt.legend()
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.C_norm, label='Cytokines')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        if ids is None:
            plt.title(f"Évolution des concentrations maximales")
        else :
            plt.title(f"Évolution des concentrations maximales. Iterations {ids}")
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.Coll_norm, label='Collagène')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        if ids is None:
            plt.title(f"Évolution des concentrations maximales")
        else :
            plt.title(f"Évolution des concentrations maximales. Iterations {ids}")
        plt.legend()
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.C_norm, label='Cytokines')
        plt.plot(self.t, self.F_norm, label='Fibroblastes')
        plt.plot(self.t, self.Coll_norm, label='Collagène')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        if ids is None:
            plt.title(f"Évolution des concentrations maximales")
        else :
            plt.title(f"Évolution des concentrations maximales. Iterations {ids}")
        plt.legend()
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.t, self.C_norm, label='Cytokines')
        plt.plot(self.t, self.F_norm, label='Fibroblastes')
        plt.plot(self.t, self.Coll_norm, label='Collagène')
        plt.yscale('log')
        plt.xlabel("Temps (jours)")
        plt.ylabel("Valeur maximale")
        if ids is None:
            plt.title(f"Évolution des concentrations maximales")
        else :
            plt.title(f"Évolution des concentrations maximales. Iterations {ids}")
        plt.legend()
        plt.legend()
        plt.grid(True)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        plt.show()
        
        
    def demie_vie_cyto(self):
        plt.plot(self.t, self.C_norm)
        plt.yscale('log')
        #plt.xscale('log')
        plt.title("Evolution du taux maximal de cytokines au cours du temps")
        plt.xlabel("temps (jour)")
        plt.ylabel("Max. cyt. (echelle log.)")
        plt.show()
    '''
    def track_chemo(self):
        if len(self.C_results) != len(self.F_results):
            raise ValueError("C_norm et F_norm doivent avoir la même longeur !")
        
        from matplotlib import cm
            
        traj_x = []
        traj_y = []
        traj_z = []
        for i in range(len(self.C_results)):
            max_mat = np.maximum(self.C_results[i], self.F_results[i])
            pos_max = np.unravel_index(np.argmax(max_mat), max_mat.shape)
            traj_x.append(pos_max[0])
            traj_y.append(pos_max[1])
            traj_z.append(pos_max[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Création d'un colormap basé sur le temps
        colors = cm.viridis(np.linspace(0, 1, len(self.C_results)))

        for i in range(len(self.C_results) - 1):
            ax.plot(traj_x[i:i+2], traj_y[i:i+2], traj_z[i:i+2], color=colors[i])

        # Scatter pour les points avec couleur selon le temps
        p = ax.scatter(traj_x, traj_y, traj_z, c=np.arange(len(self.C_results)), cmap='viridis')

        fig.colorbar(p, ax=ax, label='Temps')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Trajectoire 3D du point max colorée selon le temps')

        plt.show()
    '''
    
    def track_chemo(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import cm
    
        if len(self.C_results) != len(self.F_results):
            raise ValueError("C_norm et F_norm doivent avoir la même longueur !")
    
        traj_x = []
        traj_y = []
        traj_z = []
    
        traj_C = []
        traj_F = []
    
        for i in range(len(self.C_results)):
            # Max indépendant
            pos_C = np.unravel_index(np.argmax(self.C_results[i]), self.C_results[i].shape)
            pos_F = np.unravel_index(np.argmax(self.F_results[i]), self.F_results[i].shape)
            traj_C.append(pos_C)
            traj_F.append(pos_F)
    
            # Max commun (élément le plus élevé entre C et F)
            max_mat = np.maximum(self.C_results[i], self.F_results[i])
            pos_max = np.unravel_index(np.argmax(max_mat), max_mat.shape)
            traj_x.append(pos_max[0])
            traj_y.append(pos_max[1])
            traj_z.append(pos_max[2])
    
        # Décomposition pour tracer
        traj_C_x, traj_C_y, traj_C_z = zip(*traj_C)
        traj_F_x, traj_F_y, traj_F_z = zip(*traj_F)
    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
        # Couleurs fixes pour les deux trajectoires individuelles
        ax.plot(traj_C_x, traj_C_y, traj_C_z, color='blue', label='Max C')
        ax.plot(traj_F_x, traj_F_y, traj_F_z, color='red', label='Max F')
    
        # Colormap pour les points du max commun
        scatter = ax.scatter(traj_x, traj_y, traj_z, c=np.arange(len(self.C_results)), cmap='viridis', label='Max commun (C vs F)')
        ax.plot(traj_x, traj_y, traj_z, color='green', linestyle='--')  # ligne de liaison
    
        fig.colorbar(scatter, ax=ax, label='Temps')
    
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Trajectoires 3D')
        ax.legend()
    
        plt.show()


        

#%%

def norm(u,v):
    res = 0
    w = u - v
    for i in range(len(u)):
        res += w[i]**2
    return np.sqrt(res)

def compute_cytokines(vtk_path, mesh_path):
    mesh = Mesh(mesh_path)
    exp = Export(vtk_path)
    N = 11
    #norms = exp.get_von_mises_stresses(exp.get_tensors(name="nodal_strain"))
    norms = exp.compute_tensor_norms(exp.get_tensors(name="nodal_strain"))
    nodes = mesh.get_nodes()
    
    # Convertir en listes pour faciliter l'indexation
    node_coords = list(nodes.values())
    node_indices = list(nodes.keys())
    
    Cyt = np.zeros((N,N,N))
    x = np.linspace(0,1,N)
    y = np.linspace(0,1,N)
    z = np.linspace(0,1,N)
    
    #print(f"Nombre de nœuds: {len(node_coords)}")
    #print(f"Nombre de normes: {len(norms)}")
    
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                grid_point = np.array([x[i], y[j], z[k]])
                
                # Option 1: Trouver le nœud le plus proche
                min_distance = float('inf')
                closest_node_idx = -1
                
                for node_idx, node_coord in enumerate(node_coords):
                    node_pos = np.array([node_coord[0], node_coord[1], node_coord[2]])
                    distance = np.linalg.norm(node_pos - grid_point)
                    
                    if distance < min_distance:
                        min_distance = distance
                        closest_node_idx = node_idx
                
                # Utiliser la norme du nœud le plus proche (dans la limite de distance)
                if min_distance <= 1.0:  # ou un seuil plus petit
                    Cyt[i,j,k] = norms[closest_node_idx]
                else:
                    Cyt[i,j,k] = 0  # ou une valeur par défaut
                    
    return Cyt

#%%
if __name__ == "__main__":
    from pathlib import Path

    vtk_path = Path("Z:/Automatisation/DoNotTouch_Files/init.vtk")
    feb_path = Path("Z:/Automatisation/DoNotTouch_Files/init.feb")
    mesh_path = Path("Z:/Model/cube.msh")
    #♥fib_path = Path("Z:/Boucle_Complete/fib_final_4.txt")
    #coll_path = Path("Z:/Boucle_Complete/coll_final_4.txt")
    
    exp = Export(vtk_path)

    Cyt = compute_cytokines(vtk_path, mesh_path)
    #Fib = compute(fib_path)
    
    #Coll = compute(coll_path)


         
    simulation = Liver_Fibrosis_Solver(
        N=11,
        L=1.0,
        dt=0.01,
        T=730,
        C=Cyt,
        F=None,
        Coll=None,
        vtk_filename=vtk_path,
        feb_filename=feb_path
    )
    
    norms = simulation.load_tensors_norms()
    
    simulation.initialize()
    #simulation.plot()
    simulation.solve(store_full_results=True)
    
    # Afficher les résultats toutes les 10 itérations
    for i in range(0, simulation.n_steps, 40000):
        #if i > 0:  # Éviter l'état initial qui est déjà affiché
        simulation.plot_slice_custom(step=i, x_slice=5, y_slice=7, z_slice=10)
            
    simulation.plot_norms()
    simulation.track_chemo()
    
    save_to_file('Z:/Automatisation/results/coll_final_1.txt', simulation.x, simulation.y, simulation.z, simulation.Coll_results[-1])
    save_to_file('Z:/Automatisation/results/fib_final_1.txt', simulation.x, simulation.y, simulation.z, simulation.F_results[-1])
    
    plt.plot(simulation.t, simulation.C_norm)
    plt.yscale('log')
    #plt.xscale('log')
    plt.title("Evolution du taux maximal de cytokines au cours du temps")
    plt.xlabel("temps (jour)")
    plt.ylabel("Max. cyt. (echelle log.)")
    plt.show()
    
    simulation.demie_vie_cyto()
    
