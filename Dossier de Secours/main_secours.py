import numpy as np
import matplotlib.pyplot
from mesh import Mesh

from liver_fibrosis_solver import *
from coll_interpolate import *
from create_feb import *
from copy_feb_files import copy
from TRUE_TRUE_TRUE_5 import *


iterations = np.arange(1,5)


def run_and_save_bio_solver(solver: Liver_Fibrosis_Solver):
    solver.initialize()
    
    
if __name__ == '__main__':
    
    mesh_path = Path("Z:/Automatisation/DoNotTouch_Files/cube.msh")
    
    virgin_path = Path('Z:/Automatisation/DoNotTouch_Files/virgin.feb')
    destination_folder = Path('Z:/Automatisation/FEB_Files')
    
    mesh = Mesh(mesh_path)
    Elements = mesh.get_elements()
    Nodes = mesh.get_nodes()
    
    
    
    print('Début boucle')
    print('\n')
    for ids in iterations:
        
        
        coll_final_path = Path(f'Z:/Automatisation/results/coll_final_{ids}.txt')
        fib_final_path = Path(f'Z:/Automatisation/results/fib_final_{ids}.txt')
        coll_value_interpolate_path = Path(f"Z:/Automatisation/results/coll_value_interpolate{ids}.txt")
        feb_destination_path = Path(f'Z:/Automatisation/FEB_Files/iter{ids}.feb')
        
        if ids == 1:
            vtk_path = Path("Z:/Automatisation/DoNotTouch_Files/init.vtk")
            feb_path = Path("Z:/Automatisation/DoNotTouch_Files/init.feb")
            exp = Export(vtk_path)
            Cyt = compute_cytokines(vtk_path, mesh_path)
            
            solve = Liver_Fibrosis_Solver(
                N=11,
                L=1.0,
                dt=0.001,
                T=10,
                C=Cyt,
                F=None,
                Coll=None,
                vtk_filename=vtk_path,
                feb_filename=feb_path
            )
            
            print(iterations)
            
            norms = solve.load_tensors_norms()
            
            solve.initialize()
            solve.plot()
            solve.solve(store_full_results=True)
            
            for i in range(0, solve.n_steps, 1000):
                if i > 0:  
                    solve.plot_slice_custom(step=i, x_slice=5, y_slice=5, z_slice=10)
                    
            solve.plot_norms()
            
            save_to_file(coll_final_path, solve.x, solve.y, solve.z, solve.Coll_results[-1])
            save_to_file(fib_final_path, solve.x, solve.y, solve.z, solve.F_results[-1])
            
            ''' Fin Simulation Bio '''
            ''' Interpolation des Résultats'''
            print('Interpolation des Résultats...')
            print('\n')

            Coll = load_from_file(coll_final_path)  
            Val = build_interpolated_element_values(Elements, Nodes, Coll)
            save_coll_elmt_interpolate(coll_value_interpolate_path, Val)
            
            ''' Creation du nouveau fichier feb'''
            print('Création du nouveau fichier .feb')
            print('\n')
            
            destination_filename = f'iter{ids}.feb'
            copy(virgin_path, destination_folder, destination_filename)
            
            
            '''Récupération des collagènes interpolées sous la forme d'un dictionnaire'''
            coll_path = coll_final_path
            coll = data_dict(coll_value_interpolate_path)
            for key, value in coll.items():
                '''
                if value < 0.050 :
                    coll[key] = value * 100 + 3
                else :
                    coll[key] = 80
                '''
                coll[key] = value * 100 + 3
            print(coll)
            print('\n')
            
            '''Liste de module de Young'''
            young_list = young_liste(coll)
            print(young_list)
            print('\n')
            
            '''Initialisation d'une liste de dictionnaire'''
            dico_list = [{} for _ in range(len(young_list)-1)]
            print(dico_list)
            print('\n')
            
            for i in range(len(young_list)-1):
                a = young_list[i]
                b = young_list[i+1]
                mean = (a + b) / 2
                dico_list[i][mean] = []
                
            print(dico_list)
            
            # Remplissage
            for key, value in coll.items():
                for i in range(len(young_list)-1):
                    a = young_list[i]
                    b = young_list[i+1]
                    if a <= value < b:
                        mean = (a + b) / 2
                        dico_list[i][mean].append(key)
                        if i == 4:
                            print(mean, value)
                        break
            print('Interval : \n',dico_list,'\n')
            max_value = max(coll.values())
            if max_value >= young_list[-1]:
                mean = (young_list[-1] + (young_list[-1] + 1)) / 2  # ou tout autre borne supérieure
                dico_list.append({mean: []})
                for key, value in coll.items():
                    if value >= young_list[-1]:
                        dico_list[-1][mean].append(key)
            print('Interval : \n',dico_list,'\n')
            
            ''''Vérification du nombre d'éléments'''
            for i in range(len(dico_list)):
                for values in dico_list[i].values():
                    print('\n')
                    print(values)
                    print(len(values))
                    print('\n')
                    
        #%%      
            elmt_data = elmt_dta(dico_list, mesh_path)
            
            print('Interval : \n',dico_list,'\n')
            for i in range(len(dico_list)):
                if len(list(dico_list[i].values())[0]) > 0:
                    mat = Material(ids=f'{i+1}', name=f'Material{i+1}', E=list(dico_list[i].keys())[0])
                    insert_material_block_in(feb_destination_path, mat)
                    
                    elmt = Element(elmt_data=elmt_data[i], name=f'Part{i+1}', part=f'Part{i+1}', mat=f'{i+1}')
                    insert_element_block_in(feb_destination_path, elmt)
                    
                    solid_domain = SolidDomain(name=f'Part{i+1}', mat=f'Material{i+1}')
                    insert_solid_domain_in(feb_destination_path, solid_domain)
                
            print('\n')
            print(dico_list)
            
            ''' Fin du montage .feb'''
            ''' Run du .feb avec FEBio'''
            import subprocess
            import os
            import sys
            import time
            import pyautogui
            
            vtk_outpout_name = f'iter{ids}.vtk'
        
            powershell_cmd = f'& "C:\\Program Files\\FEBioStudio\\bin\\FEBioStudio2.exe" "{feb_destination_path.as_posix()}"'
        
            subprocess.run(["powershell", "-Command", powershell_cmd])
            
            time.sleep(5)
            
            pyautogui.press('enter')
            time.sleep(2)
            
            pyautogui.press('F5')
            time.sleep(5)
            
            pyautogui.press('enter')
            time.sleep(5)
            #pyautogui.press('left')
            #time.sleep(0.5)
            pyautogui.press('enter')
            time.sleep(5)
            pyautogui.press('enter')
            time.sleep(5)
            pyautogui.press('right')
            time.sleep(0.3)
            
            # Ctrl+Shift+S pour "Enregistrer sous"
            pyautogui.hotkey('ctrl', 'shift', 's')
            time.sleep(1)
            
            
            pyautogui.write(str(vtk_outpout_name))
            time.sleep(0.5)
            # Tab pour se positionner sur la liste des formats
            pyautogui.press('tab')
            
            time.sleep(0.5) 
    
            # 8 fois flèche bas pour aller sur VTK
            for _ in range(8):
                pyautogui.press('down')
                time.sleep(0.1)
    
            # 4 fois entrée pour valider
            for _ in range(4):
                pyautogui.press('enter')
                time.sleep(0.3)
        else:
            idd=ids-1
            vtk_path = Path(f'Z:/Automatisation/FEB_Files/jobs/iter{idd}.vtk')
            vtk_copy_path = Path(f'Z:/Automatisation/FEB_Files')
            vtk_name = f'iter{idd}.vtk'
            copy(vtk_path, vtk_copy_path, vtk_name)
            
            vtk_path = Path(f"Z:/Automatisation/FEB_Files/iter{idd}.vtk")
            feb_path = Path(f"Z:/Automatisation/FEB_Files/iter{idd}.feb")
            exp = Export(vtk_path)
            
            fib_path = Path(f'Z:/Automatisation/results/fib_final_{idd}.txt')
            coll_path = Path(f'Z:/Automatisation/results/coll_final_{idd}.txt')
            
            Cyt = compute_cytokines(vtk_path, mesh_path)
            Fib = compute(fib_path)
            Coll = compute(coll_path)
            
            solve = Liver_Fibrosis_Solver(
                N=11,
                L=1.0,
                dt=0.001,
                T=3,
                C=Cyt,
                F=Fib,
                Coll=Coll,
                vtk_filename=vtk_path,
                feb_filename=feb_path
            )
            
            print(iterations)
            
            norms = solve.load_tensors_norms()
            
            solve.initialize()
            solve.plot()
            solve.solve(store_full_results=True)
            
            for i in range(0, solve.n_steps, 1000):
                if i > 0:  
                    solve.plot_slice_custom(step=i, x_slice=5, y_slice=5, z_slice=10)
                    
            solve.plot_norms()
            
            
            save_to_file(coll_final_path, solve.x, solve.y, solve.z, solve.Coll_results[-1])
            save_to_file(fib_final_path, solve.x, solve.y, solve.z, solve.F_results[-1])
            
            ''' Fin Simulation Bio '''
            
            
            
            
            
            ''' Interpolation des Résultats'''
            print('Interpolation des Résultats...')
            print('\n')

            Coll = load_from_file(coll_final_path)  
            Val = build_interpolated_element_values(Elements, Nodes, Coll)
            save_coll_elmt_interpolate(coll_value_interpolate_path, Val)
            
            ''' Creation du nouveau fichier feb'''
            print('Création du nouveau fichier .feb')
            print('\n')
            
            destination_filename = f'iter{ids}.feb'
            copy(virgin_path, destination_folder, destination_filename)
            
            
            '''Récupération des collagènes interpolées sous la forme d'un dictionnaire'''
            coll_path = coll_final_path
            coll = data_dict(coll_value_interpolate_path)
            for key, value in coll.items():
                '''
                if value < 0.050 :
                    coll[key] = value * 100 + 3
                else :
                    coll[key] = 80
                '''
                coll[key] = value * 100 + 3
            print(coll)
            print('\n')
            
            '''Liste de module de Young'''
            young_list = young_liste(coll)
            print(young_liste)
            print('\n')
            
            '''Initialisation d'une liste de dictionnaire'''
            dico_list = [{} for _ in range(len(young_list)-1)]
            print(dico_list)
            print('\n')
            
            for i in range(len(young_list)-1):
                a = young_list[i]
                b = young_list[i+1]
                mean = (a + b) / 2
                dico_list[i][mean] = []
                
            print(dico_list)
            
            # Remplissage
            for key, value in coll.items():
                for i in range(len(young_list)-1):
                    a = young_list[i]
                    b = young_list[i+1]
                    if a <= value < b:
                        mean = (a + b) / 2
                        dico_list[i][mean].append(key)
                        if i == 4:
                            print(mean, value)
                        break
            print('Interval : \n',dico_list,'\n')
            max_value = max(coll.values())
            if max_value >= young_list[-1]:
                mean = (young_list[-1] + (young_list[-1] + 1)) / 2  # ou tout autre borne supérieure
                dico_list.append({mean: []})
                for key, value in coll.items():
                    if value >= young_list[-1]:
                        dico_list[-1][mean].append(key)
            print('Interval : \n',dico_list,'\n')
            
            ''''Vérification du nombre d'éléments'''
            for i in range(len(dico_list)):
                for values in dico_list[i].values():
                    print('\n')
                    print(values)
                    print(len(values))
                    print('\n')
                    
        #%%      
            elmt_data = elmt_dta(dico_list, mesh_path)
            
            print('Interval : \n',dico_list,'\n')
            for i in range(len(dico_list)):
                if len(list(dico_list[i].values())[0]) > 0:
                    mat = Material(ids=f'{i+1}', name=f'Material{i+1}', E=list(dico_list[i].keys())[0])
                    insert_material_block_in(feb_destination_path, mat)
                    
                    elmt = Element(elmt_data=elmt_data[i], name=f'Part{i+1}', part=f'Part{i+1}', mat=f'{i+1}')
                    insert_element_block_in(feb_destination_path, elmt)
                    
                    solid_domain = SolidDomain(name=f'Part{i+1}', mat=f'Material{i+1}')
                    insert_solid_domain_in(feb_destination_path, solid_domain)
                
            print('\n')
            print(dico_list)
            
            
            ''' Fin du montage .feb'''
            
            
            ''' Run du .feb avec FEBio'''
            import subprocess
            import os
            import sys
            import time
            import pyautogui
            
            vtk_outpout_name = f'iter{ids}.vtk'

            powershell_cmd = f'& "C:\\Program Files\\FEBioStudio\\bin\\FEBioStudio2.exe" "{feb_destination_path.as_posix()}"'
        
            subprocess.run(["powershell", "-Command", powershell_cmd])
            
            time.sleep(5)
            
            pyautogui.press('enter')
            time.sleep(2)
            
            pyautogui.press('F5')
            time.sleep(5)
            
            pyautogui.press('enter')
            time.sleep(5)
            #pyautogui.press('left')
            #time.sleep(0.5)
            pyautogui.press('enter')
            time.sleep(5)
            pyautogui.press('enter')
            time.sleep(5)
            pyautogui.press('right')
            time.sleep(0.3)
            
            # Ctrl+Shift+S pour "Enregistrer sous"
            pyautogui.hotkey('ctrl', 'shift', 's')
            time.sleep(1)
            
            pyautogui.write(str(vtk_outpout_name))
            time.sleep(0.5)
            # Tab pour se positionner sur la liste des formats
            pyautogui.press('tab')
            
            time.sleep(0.5) 
    
            # 8 fois flèche bas pour aller sur VTK
            for _ in range(8):
                pyautogui.press('down')
                time.sleep(0.1)
    
            # 4 fois entrée pour valider
            for _ in range(4):
                pyautogui.press('enter')
                time.sleep(0.3)
    
    print('Fin boucle')
    print('\n')
            