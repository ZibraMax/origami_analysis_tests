from typing import List
from .Elements import Panel, Bar, Hinge, RectangularPanel, TriangularPanel
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


class Opening(Bar):
    def __init__(self, nodes):
        super().__init__(nodes)

    def set_domanin(self, domain):
        super().set_domanin(domain)

    def is_interior(self, point, tol=1e-8):
        A = self.coords[0]
        B = self.coords[1]

        point = np.array(point)
        AB = B - A
        AP = point - A
        if np.linalg.norm(AP) < tol or np.linalg.norm(B - point) < tol:
            return False  # Returns false if point is on the edge

        if np.allclose(AB, 0, atol=tol):
            return np.allclose(AP, 0, atol=tol)
        cross = np.cross(AB, AP)
        return np.linalg.norm(cross) < tol


class Geometry():
    def __init__(self, ngdl_per_node=6):
        self.base_nodes = np.zeros((0, 3), dtype=float)
        self.base_gdls = np.zeros((0, 6), dtype=int)
        self.ngdl_per_node: int = ngdl_per_node
        self.panels: List[Panel] = []
        self.bars: List[Bar] = []
        self.hinges: List[Hinge] = []
        self.openings: List[Opening] = []
        self.meshed = False
        self.ebc = []
        self.nbc = []
        self.tie_nodes = []

        self.nodes = np.zeros((0, 3), dtype=float)
        self.gdls = np.zeros((0, 6), dtype=int)
        self.node_map = {}
        self.dictionary = []
        self.types = []
        self.properties = {}

    def test_point_vertices(self, vertex, tol=1e-10):
        if len(self.nodes) == 0:
            return False, None
        if not isinstance(vertex, np.ndarray):
            vertex = np.array(vertex)
        distances = np.linalg.norm(self.nodes-vertex, axis=1)
        if np.min(distances) < tol:
            return True, int(np.argmin(distances))
        return False, None

    def test_points_vertices(self, point, tol=1e-10):
        if len(self.nodes) == 0:
            return []
        if not isinstance(point, np.ndarray):
            point = np.array(point)
        indices = []
        for i, v in enumerate(self.nodes):
            if np.linalg.norm(v-point) < tol:
                indices.append(i)
        return indices

    def get_index_duplicates(self, tol=1e-10):
        if len(self.nodes) == 0:
            return {}
        duplicates = {}
        for i, v in enumerate(self.nodes):
            indices = self.test_points_vertices(v, tol=tol)
            if len(indices) > 1:
                duplicates[i] = indices
        return duplicates

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        return ax

    def add_node(self, coord):
        n = len(self.base_gdls)*self.ngdl_per_node
        gdl = np.arange(n, n+self.ngdl_per_node)
        self.base_nodes = np.vstack([self.base_nodes, coord])
        self.base_gdls = np.vstack([self.base_gdls, gdl])

    def add_nodes(self, coords):
        for coord in coords:
            self.add_node(coord)

    def add_panel(self, panel: Panel):
        panel.set_domanin(self)
        self.panels.append(panel)

    def add_bar(self, bar: Bar):
        bar.set_domanin(self)
        self.bars.append(bar)

    def add_hinge(self, hinge: Hinge):
        hinge.set_domanin(self)
        self.hinges.append(hinge)

    def add_opening(self, opening: Opening):
        opening.set_domanin(self)
        self.openings.append(opening)

    def add_ebc(self, node, values=None):
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before adding boundary conditions.")
        if values is None:
            values = [1]*self.ngdl_per_node
        if len(values) != self.ngdl_per_node:
            values = values + [0]*(self.ngdl_per_node - len(values))
        if self.node_map.get(node) is None:
            raise ValueError(f"Node {node} not found in geometry.")
        self.ebc.append([self.node_map[node][0], *values])

    def add_nbc(self, node, values=None):
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before adding boundary conditions.")
        if values is None:
            values = [0]*self.ngdl_per_node
        if len(values) != self.ngdl_per_node:
            values = values + [0]*(self.ngdl_per_node - len(values))
        if self.node_map.get(node) is None:
            raise ValueError(f"Node {node} not found in geometry.")
        self.nbc.append([self.node_map[node][0], *values])

    def add_bc_plane(self, plane='z', coord=0.0, values=None, tol=1e-8):
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before adding boundary conditions.")
        if plane not in ['x', 'y', 'z']:
            raise ValueError("Plane must be 'x', 'y', or 'z'.")
        axis = {'x': 0, 'y': 1, 'z': 2}[plane]
        for i, node in enumerate(self.nodes):
            if abs(node[axis]-coord) < tol:
                self.ebc.append([i, *values])

    def add_load_plane(self, plane='z', coord=0.0, values=None, tol=1e-8):
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before adding boundary conditions.")
        if plane not in ['x', 'y', 'z']:
            raise ValueError("Plane must be 'x', 'y', or 'z'.")
        axis = {'x': 0, 'y': 1, 'z': 2}[plane]
        for i, node in enumerate(self.base_nodes):
            if abs(node[axis]-coord) < tol:
                self.nbc.append([self.node_map[i][0], *values])

    def get_nodes_plane(self, plane='z', coord=0.0, tol=1e-8, basenodes=False):
        lookupnodes = self.base_nodes if basenodes else self.nodes
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before adding boundary conditions.")
        if plane not in ['x', 'y', 'z']:
            raise ValueError("Plane must be 'x', 'y', or 'z'.")
        axis = {'x': 0, 'y': 1, 'z': 2}[plane]
        nodes = []
        for i, node in enumerate(lookupnodes):
            if abs(node[axis]-coord) < tol:
                nodes.append(i)
        return nodes

    def add_tie_nodes(self, node1, node2):
        if dofs is None:
            dofs = list(range(self.ngdl_per_node))
        self.tie_nodes.append((node1, node2))

    def mesh(self, n=1, mesh_hinges=True):
        if self.meshed:
            return

        self.nodes = []
        elements = []
        tie_nodes = []
        types = []
        properties = {}
        properties["th"] = []

        for panel in self.panels:
            for hinge in self.hinges:
                i1, i2, i3, i4 = hinge.nodes
                if i2 in panel.nodes and i3 in panel.nodes:
                    hinge.panels.append(panel)

            res = panel.mesh(n=n)
            sub_nodes = panel.discretized_nodes
            sub_elements = panel.discretized_elements
            sub_type = "Shell"
            sub_thickness = panel.thickness
            index_offset = len(self.nodes)
            for node in sub_nodes:
                self.nodes.append(node)
            for elem in sub_elements:
                a = 0
                new_elem = [idx + index_offset for idx in elem]
                elem[:] = [int(i) for i in new_elem]
                elements.append(elem)
                types.append(sub_type)
                properties["th"].append(sub_thickness)

        self.nodes = np.array(self.nodes)

        for i, hinge in enumerate(self.hinges):
            hinge.mesh(n)
            for subelement in hinge.discretized_elements:
                idxn1, idxn2, idxn3, idxn4 = subelement
                if None in [idxn1, idxn2, idxn3, idxn4]:
                    raise ValueError(
                        f"Hinge nodes not found in new nodes list: {hinge.nodes}")
                elements.append([idxn1, idxn2, idxn3, idxn4])
                # hinge.discretized_nodes = [idxn1, idxn2, idxn3, idxn4]
                types.append("OriHinge")

        for bar in self.bars:
            n1, n2 = bar.nodes
            v1 = self.base_nodes[n1]
            v2 = self.base_nodes[n2]
            idxn1 = self.test_point_vertices(v1)[1]
            idxn2 = self.test_point_vertices(v2)[1]
            if None in [idxn1, idxn2]:
                raise ValueError(
                    f"Interior bar nodes not found in new nodes list: {n1},{n2}")
            elements.append([idxn1, idxn2])
            bar.discretized_nodes = [idxn1, idxn2]
            types.append("Bar")

        duplicates = self.get_index_duplicates()
        for parent, duplicate_indices in duplicates.items():
            test_node = self.nodes[duplicate_indices[0]]
            flag = True
            for opening in self.openings:
                if opening.is_interior(test_node):
                    flag = False
                    break
            if flag:
                for idx in duplicate_indices[1:]:
                    tie_nodes.append((duplicate_indices[0], idx))
        self.dictionary = elements
        self.tie_nodes = tie_nodes
        self.types = types
        properties["problem"] = "ShellAndHinge"
        self.properties = properties
        self.gdls = np.zeros((len(self.nodes), self.ngdl_per_node), dtype=int)

        self.node_map = {}
        for i, node in enumerate(self.base_nodes):
            idxs = self.test_points_vertices(node)
            self.node_map[i] = idxs

        self.meshed = True

    @staticmethod
    def from_json_file(file_path, t=0.2):
        data = json.load(open(file_path, 'r'))
        return Geometry.from_json(data, t=t)

    @staticmethod
    def from_json(data, t=0.2):
        O = Geometry(ngdl_per_node=6)
        nodes = np.array(data['nodes'])
        dictionary = data['dictionary']
        types = data['types']
        for node in nodes:
            O.add_node(node)
        for elem_nodes, elem_type in zip(dictionary, types):
            if elem_type == 'Shell' or "ASDShell" in elem_type:
                th = data["properties"].get("thickness", t)
                if len(elem_nodes) == 3:
                    panel = TriangularPanel(elem_nodes, thickness=th)
                elif len(elem_nodes) == 4:
                    panel = RectangularPanel(elem_nodes, thickness=th)
                O.add_panel(panel)
            elif elem_type == 'Bar' or elem_type == 'L1V':
                bar = Bar(elem_nodes)
                O.add_bar(bar)
            elif elem_type == 'OriHinge':
                t1 = data["properties"].get("theta1", R(10))
                t2 = data["properties"].get("theta2", R(350))
                hinge = Hinge(elem_nodes, theta1=t1, theta2=t2)
                O.add_hinge(hinge)
            elif elem_type == 'Opening':
                opening = Opening(elem_nodes)
                O.add_opening(opening)
            else:
                raise ValueError(f"Unknown element type: {elem_type}")
        return O

    def to_json(self):
        dicttypes = {"Shell": "T1V", "Bar": "L1V",
                     "OriHinge": "OH", "Opening": "OH"}
        if not self.meshed:
            raise ValueError(
                "Geometry must be meshed before exporting to JSON.")
        data = {}
        data['nodes'] = self.nodes.tolist()
        data['dictionary'] = self.dictionary
        data['types'] = [dicttypes[t] for t in self.types]
        for i in range(len(data["dictionary"])):
            if len(data["dictionary"][i]) == 3 and data["types"][i] == "T1V":
                data["types"][i] = "T1V"
            elif len(data["dictionary"][i]) == 4 and data["types"][i] == "T1V":
                data["types"][i] = "C1V"
        data['properties'] = self.properties
        return data
