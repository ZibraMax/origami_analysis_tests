import numpy as np
from .Elements import Panel, Bar, Hinge


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
        if np.allclose(AB, 0, atol=tol):
            return np.allclose(AP, 0, atol=tol)
        cross = np.cross(AB, AP)
        return np.linalg.norm(cross) < tol


class Geometry():
    def __init__(self, ngdl_per_node=6):
        self.base_nodes = np.zeros((0, 3), dtype=float)
        self.base_gdls = np.zeros((0, 6), dtype=int)
        self.ngdl_per_node: int = ngdl_per_node
        self.panels: Panel = []
        self.bars: Bar = []
        self.hinges: Hinge = []
        self.openings: Opening = []

        self.nodes = np.zeros((0, 3), dtype=float)
        self.gdls = np.zeros((0, 6), dtype=int)

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
        if not isinstance(vertex, np.ndarray):
            vertex = np.array(vertex)
        indices = []
        for i, v in enumerate(self.nodes):
            if np.linalg.norm(v-point) < tol:
                indices.append(i)
        return indices

    def get_index_duplicates(self, tol=1e-10):
        if len(self.nodes) == 0:
            return []
        duplicates = {}
        for i, v in enumerate(self.nodes):
            indices = self.test_points_vertices(v, tol=tol)
            if len(indices) > 1:
                duplicates[i] = indices
        return duplicates

    def add_node(self, coord):
        n = len(self.base_gdls)*self.ngdl_per_node
        gdl = np.arange(n, n+self.ngdl_per_node)
        self.base_nodes = np.vstack([self.base_nodes, coord])
        self.base_gdls = np.vstack([self.base_gdls, gdl])

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
