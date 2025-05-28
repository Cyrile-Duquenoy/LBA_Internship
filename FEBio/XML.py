import numpy as np
from mesh import Mesh
from pathlib import Path
from xml.etree.ElementTree import Element, SubElement, ElementTree, parse
import xml.etree.ElementTree as ET

from classes import (Material,
                     Element,
                     Node,
                     SolidDomain)

def create_material_block(material: Material):
    """Crée un élément <material> XML complet pour FEBio."""
    mat_block = ET.Element("material", {
        "id": str(material.ids),
        "type": material.types,
        "name": material.name 
    })

    # Ajouter les propriétés internes
    e_tag = ET.SubElement(mat_block, "E")
    e_tag.text = str(material.E)

    v_tag = ET.SubElement(mat_block, "v")
    v_tag.text = str(material.v) if hasattr(material, "v") else "0.3"  # valeur par défaut

    d_tag = ET.SubElement(mat_block, "density")
    d_tag.text = str(material.density) if hasattr(material, "density") else "1"  # valeur par défaut

    return mat_block

def create_element_block(element: Element):
    elem_block = ET.Element("Elements", {
        'type': element.type,
        'name': element.name,
        'part': element.part,
        'mat': element.mat
        })
    return elem_block
    
    
#%%

 # Indentation propre
def indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for child in elem:
            indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = i

#%%
def insert_material_block_in(feb_path: Path, material: Material):
    """
    Insère un bloc <Elements> dans la section <Mesh> d’un fichier FEBio.
    """
    tree = parse(feb_path)
    root = tree.getroot()

    # Trouver ou créer la section Mesh
    mesh_section = root.find("Material")
    if mesh_section is None:
        mesh_section = SubElement(root, "Material")

    # Créer le bloc et l’ajouter
    block = create_material_block(material)
    mesh_section.append(block)

    indent(root)
    tree.write(feb_path, encoding="utf-8", xml_declaration=True)
    print(f"✅ Bloc <Elements> ajouté dans : {feb_path}")
    
    
#%%

def insert_element_block_in(feb_path: Path, element):
    """
    Insère un élément (volume hex8) dans un fichier .feb FEBio, avec indentation propre.
    """
    tree = parse(feb_path)
    root = tree.getroot()

    # Essayer de trouver ou créer la section Geometry
    geometry_section = root.find("Mesh")
    if geometry_section is None:
        geometry_section = root.find(".//Mesh")
    if geometry_section is None:
        geometry_section = SubElement(root, "Mesh")

    # Créer la balise Elements
    elements_section = create_element_block(element)

    # Ajout de chaque élément avec ID
    for elmt_id, node_list in element.elmt_data.items():
        node_str = ",".join(map(str, node_list))
        elem_tag = SubElement(elements_section, "elem", {"id": str(elmt_id)})
        elem_tag.text = node_str
    geometry_section.append(elements_section)

    indent(root)

    # Sauvegarde propre
    tree.write(feb_path, encoding='utf-8', xml_declaration=True)
    print(f"✅ Élément '{element.name}' inséré dans : {feb_path}")
    
#%%

def create_node_block(node:Node):
    node_block = ET.Element("Nodes", {
        'name': node.name
        })
    return node_block

def insert_node_block_in(feb_path: Path, node):
    tree = parse(feb_path)
    root = tree.getroot()
    
    geometry_section = root.find("Mesh")
    if geometry_section is None:
        geometry_section = root.find(".//Mesh")
    if geometry_section is None:
        geometry_section = SubElement(root, "Mesh")
    
    # Créer la balise Nodes
    node_section = create_node_block(node)

    # Ajout de chaque élément avec ID
    for node_id, coord_list in node.elmt_data.items():
        coord_str = ",".join(map(str, coord_list))
        elem_tag = SubElement(node_section, "elem", {"id": str(elmt_id)})
        elem_tag.text = coord_str
    geometry_section.append(node_section)

    indent(root)

    # Sauvegarde propre
    tree.write(feb_path, encoding='utf-8', xml_declaration=True)
    print(f"✅ Nodes '{node.name}' inséré dans : {feb_path}")
    
#%%

def create_solid_domain_block(solid_domain):
    solid_domain_block = ET.Element("SolidDomain", {
        'name': solid_domain.name,
        'mat': solid_domain.mat
        })
    return solid_domain_block

def insert_solid_domain_in(feb_path: Path, solid_domain):
    tree = parse(feb_path)
    root = tree.getroot()
    
    section = root.find("MeshDomains")
    if section is None:
        section = root.find(".//MeshDomains")
    if section is None:
        section = SubElement(root, "MeshDomains")
        
    solid_domain_element = create_solid_domain_block(solid_domain)
    
    section.append(solid_domain_element)
    
    tree.write(feb_path, encoding='utf-8', xml_declaration=True)

    
    
#%%
if __name__ == '__main__':
    feb_path = Path(__file__).parent / 'Model - Copie.feb'
    
    # Créer le fichier Model.feb si inexistant ou vide
    if not feb_path.exists() or feb_path.stat().st_size == 0:
        root = ET.Element("febio_spec", {"version": "4.0"})
        tree = ET.ElementTree(root)
        tree.write(feb_path, encoding="utf-8", xml_declaration=True)
        print(f"✅ Fichier initialisé : {feb_path}")
    
    mat = Material(ids=1, name='Material1', E=30)
    insert_material_block_in(feb_path, mat)
    
    elmt = Element(elmt_data={1:[1,2,3], 2:[4,5,6]}, name='Part1', part='Part1', mat='1' )
    insert_element_block_in(feb_path, elmt)
    
    
    



    
    


