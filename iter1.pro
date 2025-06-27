import numpy as np
import matplotlib.pyplot as plt
from mesh import Mesh
from liver_fibrosis_solver import *
from coll_interpolate import *
from create_feb import *
from copy_feb_files import *
from TRUE_TRUE_TRUE_5 import *
from auto_feb import *
from util import * 
            
            
if __name__ =='__main__':
    nb_iter = 4 
    nb_mat = 10 
    iteration = np.arange(1,nb_iter)
    
    mesh_path = Path('Z:/Automatisation/DoNotTouch_Files/cube.msh')
    virgin_path = Path('Z:/Automatisation/DoNotTouch_Files/virgin.feb')
    mesh_path = Path('Z:/Automatisation/DoNotTouch_Files/cube.msh')
    virgin_path = Path('Z:/Automatisation/DoNotTouch_Files/virgin.feb')
    destination_folder = Path('Z:/Automatisation/FEB_Files')
    
    mesh = Mesh(mesh_path)
    elements = mesh.get_elements()
    nodes = mesh.get_nodes()
    
    young_init_path = Path('Z:/Automatisation/DoNotTouch_Files/young_init.txt')
    
    for ids in iteration:
        #idd = ids - 1
        
        if ids == 1:
            
            young_init = read_young_list_from(young_init_path)[0]
            print(young_init)
            
            coll_final_path = Path(f'Z:/Automatisation/results/coll_final_{ids}.txt')
            fib_final_path = Path(f'Z:/Automatisation/results/fib_final_{ids}.txt')
            coll_value_interpolate_path = Path(f"Z:/Automatisation/results/coll_value_interpolate{ids}.txt")
            feb_destination_path = Path(f'Z:/Automatisation/FEB_Files/iter{ids}.feb')
            young_path = Path(f'Z:/Automatisation/results/young{ids}.txt')
            
            vtk_path = Path('Z:/Automatisation/DoNotTouch_Files/init.vtk')
            feb_path = Path('Z:/Automatisation/DoNotTouch_Files/init.feb')
            
            exp = Export(vtk_path)
            Cyt = compute_cytokines(vtk_path, mesh_path)
            solve = Liver_Fibrosis_Solver(
                N=11,
                L=1.0,
                dt=0.01,
                T=12,
                C=Cyt,
                F=None,
                Coll=None,
                vtk_filename=vtk_path,
                feb_filename=feb_path
            )
            
            exp = Export(vtk_path)
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
            
            print(' Fin Simulation Bio \n ')
            
            ''' Interpolation des Résultats'''
            print('Interpolation des Résultats...')
            print('\n')

            Coll = load_from_file(coll_final_path)  
            Val = build_interpolated_element_values(elements, nodes, Coll)
            save_coll_elmt_interpolate(coll_value_interpolate_path, Val)
            
            print(Val)
            
            print('Interpolation terminé ! \n ')
            
            ''' Creation du nouveau fichier feb'''
            print('Creation du nouveau fichier feb... \n')
            
            destination_filename = f'iter{ids}.feb'
            copy(virgin_path, destination_folder, destination_filename)
            
            print('Nouveau fichier .feb créer ! \n')
            
            new_young_dict = {}
            coll_dict = data_dict(coll_value_interpolate_path)
            for key, value in coll_dict.items():
                new_young_dict[key] = coll_dict[key] + young_init
            print(new_young_dict)
            
            y_min = min(new_young_dict.values())
            y_max = max(new_young_dict.values())
            new_young_list = np.linspace(y_min, y_max, nb_mat + 1)
            print(new_young_list)
            print(len(new_young_list))
            
            Y = []
            for i in range(len(new_young_list) - 1):
                a = new_young_list[i]
                b = new_young_list[i+1]
                mean = (a + b) / 2
                Y.append(mean)
            print('young_list : ', Y, '\n')
            
            E = []
            for key, value in new_young_dict.items():
                if new_young_list[0] <= value <= new_young_list[1]:
                    E.append(key)
            print(E)
            
            dico = [{} for _ in range(len(new_young_list) - 1)]
            for i in range(len(new_young_list) - 1):
                a = new_young_list[i]
                b = new_young_list[i+1]
                mean = (a + b) / 2
                dico[i][mean] = []
                
                for key, value in new_young_dict.items():
                    if new_young_list[i] <= value <= new_young_list[i+1]:
                        dico[i][mean].append(key)
                        new_young_dict[key] = mean
                print(len(dico[i][mean]))
            print(dico)
            
            
            young_dict_path = Path(f'Z:/Automatisation/results/young_dict{ids}.txt')
            save_young_list_to(Y, young_path)
            save_young_dict_to(new_young_dict, young_dict_path)
            
            elmt_data = elmt_dta(dico, mesh_path)
            
            print('Interval : \n',dico,'\n')
            for i in range(nb_mat):
                if len(list(dico[i].values())[0]) > 0:
                    mat = Material(ids=f'{i+1}', name=f'Material{i+1}', E=list(dico[i].keys())[0])
                    insert_material_block_in(feb_destination_path, mat)
                    
                    elmt = Element(elmt_data=elmt_data[i], name=f'Part{i+1}', part=f'Part{i+1}', mat=f'{i+1}')
                    insert_element_block_in(feb_destination_path, elmt)
                    
                    solid_domain = SolidDomain(name=f'Part{i+1}', mat=f'Material{i+1}')
                    insert_solid_domain_in(feb_destination_path, solid_domain)
                
            feb_destination_path = Path(f'Z:/Automatisation/FEB_Files/iter{ids}.feb')
            auto_feb(ids, feb_destination_path)
            
        if ids > 1:
            idd = ids - 1
            feb_destination_path = Path(f'Z:/Automatisation/FEB_Files/iter{ids}.feb')
            young_dict_path = Path(f'Z:/Automatisation/results/young_dict{idd}.txt')
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
                dt=0.01,
                T=3,
                C=Cyt,
                F=Fib,
                Coll=Coll,
                vtk_filename=vtk_path,
                feb_filename=feb_path
            )
            
            print(iteration)
            
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
            print('Fin simulation Bio ! \n')
            
            ''' Interpolation des Résultats'''
            print('Interpolation des Résultats...')
            print('\n')

            Coll = load_from_file(coll_final_path)  
            Val = build_interpolated_element_values(elements, nodes, Coll)
            save_coll_elmt_interpolate(coll_value_interpolate_path, Val)
            
            print(Val)
            
            print('Interpolation terminé ! \n ')
            
            ''' Creation du nouveau fichier feb'''
            print('Creation du nouveau fichier feb... \n')
            
            destination_filename = f'iter{ids}.feb'
            copy(virgin_path, destination_folder, destination_filename)
            
            print('Nouveau fichier .feb créer ! \n')


            new_young_dict = {}
            old_young_dict = read_young_dict_from(young_dict_path)
            print('Old young Dict : ', old_young_dict, '\n')
            coll_dict = data_dict(coll_value_interpolate_path)
            for key, value in coll_dict.items():
                new_young_dict[key] = old_young_dict[key] + coll_dict[key]
            print(new_young_dict)

            y_min = min(new_young_dict.values())
            y_max= max(new_young_dict.values())
            new_young_list = np.linspace(y_min, y_max, nb_mat + 1)
            print(new_young_list)
            print(len(new_young_list))
            
            Y = []
            for i in range(len(new_young_list) - 1):
                a = new_young_list[i]
                b = new_young_list[i+1]
                mean = (a + b) / 2
                Y.append(mean)
            print('Young_lsit : ', Y, '\n')
            
            E = []
            for key, value in new_young_dict.items():
                if new_young_list[0] <= value <= new_young_list[1]:
                    E.append(key)
            print(E)
            
            dico = [{} for _ in range(len(new_young_list) - 1)]
            for i in range(len(new_young_list) - 1):
                a = new_young_list[i]
                b = new_young_list[i+1]
                mean = (a + b) / 2
                dico[i][mean] = []
                
                for key, value in new_young_dict.items():
                    if new_young_list[i] <= value <= new_young_list[i+1]:
                        dico[i][mean].append(key)
                        new_young_dict[key] = mean
                print(len(dico[i][mean]))
            print(dico)
            
            
            young_dict_path = Path(f'Z:/Automatisation/results/young_dict{ids}.txt')
            save_young_list_to(Y, young_path)
            save_young_dict_to(new_young_dict, young_dict_path)
            
            elmt_data = elmt_dta(dico, mesh_path)
            
            print('Interval : \n',dico,'\n')
            for i in range(nb_mat):
                if len(list(dico[i].values())[0]) > 0:
                    mat = Material(ids=f'{i+1}', name=f'Material{i+1}', E=list(dico[i].keys())[0])
                    insert_material_block_in(feb_destination_path, mat)
                    
                    elmt = Element(elmt_data=elmt_data[i], name=f'Part{i+1}', part=f'Part{i+1}', mat=f'{i+1}')
                    insert_element_block_in(feb_destination_path, elmt)
                    
                    solid_domain = SolidDomain(name=f'Part{i+1}', mat=f'Material{i+1}')
                    insert_solid_domain_in(feb_destination_path, solid_domain)
                
            
            auto_feb(ids, feb_destination_path)
