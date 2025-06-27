import numpy as np
import matplotlib.pyplot as plt
import meshio


class Mesh:
    def __init__(self, filename):
        self.filename = filename
        self.nodes = {}      # node_id -> (x, y, z)
        self.elements = {}   # element_id -> [node_ids]
        self._read_msh()

    def _read_msh(self):
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line == "$Nodes":
                i += 1
                header = list(map(int, lines[i].strip().split()))
                num_entity_blocks, num_nodes, _, _ = header
                i += 1

                for _ in range(num_entity_blocks):
                    entity_header = list(map(int, lines[i].strip().split()))
                    entity_dim, entity_tag, parametric, num_nodes_in_block = entity_header
                    i += 1

                    node_ids = []
                    for _ in range(num_nodes_in_block):
                        node_ids.append(int(lines[i].strip()))
                        i += 1

                    for node_id in node_ids:
                        coords = tuple(map(float, lines[i].strip().split()))
                        self.nodes[node_id] = coords
                        i += 1

                assert lines[i].strip() == "$EndNodes"
                i += 1

            elif line == "$Elements":
                i += 1
                header = list(map(int, lines[i].strip().split()))
                num_entity_blocks, num_elements, _, _ = header
                i += 1

                for _ in range(num_entity_blocks):
                    entity_header = list(map(int, lines[i].strip().split()))
                    entity_dim, entity_tag, element_type, num_elems_in_block = entity_header
                    i += 1

                    for _ in range(num_elems_in_block):
                        parts = list(map(int, lines[i].strip().split()))
                        element_id = parts[0]
                        node_ids = parts[1:]
                        self.elements[element_id] = node_ids
                        #self.elements[element_id] = [nid - 1 for nid in node_ids]
                        i += 1

                assert lines[i].strip() == "$EndElements"
                i += 1
            else:
                i += 1

    def get_nodes(self):
        return self.nodes

    def get_elements(self):
        return self.elements



class Mesh2:
    def __init__(self, filename):
        self.filename = filename
        self.nodes = {}      # new_node_id -> (x, y, z)
        self.elements = {}   # element_id -> [new_node_ids]
        self._read_msh()

    def _read_msh(self):
        old_to_new_node_ids = {}  # old_id -> new_id

        with open(self.filename, 'r') as f:
            lines = f.readlines()

        i = 0
        node_counter = 0
        element_counter = 0

        while i < len(lines):
            line = lines[i].strip()

            if line == "$Nodes":
                i += 1
                header = list(map(int, lines[i].strip().split()))
                num_entity_blocks, num_nodes, _, _ = header
                i += 1

                for _ in range(num_entity_blocks):
                    entity_header = list(map(int, lines[i].strip().split()))
                    entity_dim, entity_tag, parametric, num_nodes_in_block = entity_header
                    i += 1

                    node_ids = []
                    for _ in range(num_nodes_in_block):
                        node_id = int(lines[i].strip())
                        node_ids.append(node_id)
                        i += 1

                    for node_id in node_ids:
                        coords = tuple(map(float, lines[i].strip().split()))
                        old_to_new_node_ids[node_id] = node_counter
                        self.nodes[node_counter] = coords
                        node_counter += 1
                        i += 1

                assert lines[i].strip() == "$EndNodes"
                i += 1

            elif line == "$Elements":
                i += 1
                header = list(map(int, lines[i].strip().split()))
                num_entity_blocks, num_elements, _, _ = header
                i += 1

                for _ in range(num_entity_blocks):
                    entity_header = list(map(int, lines[i].strip().split()))
                    entity_dim, entity_tag, element_type, num_elems_in_block = entity_header
                    i += 1

                    for _ in range(num_elems_in_block):
                        parts = list(map(int, lines[i].strip().split()))
                        element_id = element_counter
                        old_node_ids = parts[1:]
                        new_node_ids = [old_to_new_node_ids[nid] for nid in old_node_ids]
                        self.elements[element_id] = new_node_ids
                        element_counter += 1
                        i += 1

                assert lines[i].strip() == "$EndElements"
                i += 1
            else:
                i += 1
                
    def get_nodes(self):
        return self.nodes

    def get_elements(self):
        return self.elements