import numpy as np
from typing import Tuple
def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


class Element():
    def __init__(self, nodes):
        self.nodes = nodes
        self.coords = []
        self.gdl = []
        self.discretized_nodes = []

    def set_domanin(self, domain):
        self.coords = domain.base_nodes[self.nodes]
        self.gdl = domain.base_gdls[self.nodes]

    def test_point_vertices(self, vertex, nodes, tol=1e-10):
        if len(nodes) == 0:
            return False, None
        if not isinstance(vertex, np.ndarray):
            vertex = np.array(vertex)
        distances = np.linalg.norm(nodes-vertex, axis=1)
        if np.min(distances) < tol:
            return True, int(np.argmin(distances))
        return False, None


class Bar(Element):
    def __init__(self, nodes):
        super().__init__(nodes)

    def length(self):
        start, end = np.array(self.coords)
        return np.linalg.norm(end - start)


class Hinge(Element):
    def __init__(self, nodes, theta1=R(10), theta2=R(350)):
        super().__init__(nodes)
        self.theta1 = theta1
        self.theta2 = theta2
        self.panels: list[Panel] = []

    def mesh(self, n) -> Tuple[np.ndarray, list]:
        v1 = self.coords[0]
        v2 = self.coords[1]
        v3 = self.coords[2]
        v4 = self.coords[3]
        length = np.linalg.norm(v3 - v2)
        axis = (v3 - v2) / length
        origin = v2
        dt = 1 / n
        self.discretized_elements = []
        for i in range(n):
            n1 = origin + i*dt*length*axis
            n2 = origin + (i+1)*dt*length*axis
            ele_panel1 = None
            coords_panel1 = None
            idx1p1 = None
            idx2p1 = None
            ele_panel2 = None
            coords_panel2 = None
            idx1p2 = None
            idx2p2 = None
            for j, element in enumerate(self.panels[0]._discretized_elements):
                elecoords = self.panels[0].discretized_nodes[element]
                res1, idx1 = self.test_point_vertices(n1, elecoords)
                res2, idx2 = self.test_point_vertices(n2, elecoords)
                if res1 and res2:
                    ele_panel1 = self.panels[0].discretized_elements[j]
                    coords_panel1 = elecoords
                    idx1p1 = idx1
                    idx2p1 = idx2
                    break
            for j, element in enumerate(self.panels[1]._discretized_elements):
                elecoords = self.panels[1].discretized_nodes[element]
                res1, idx1 = self.test_point_vertices(n1, elecoords)
                res2, idx2 = self.test_point_vertices(n2, elecoords)
                if res1 and res2:
                    ele_panel2 = self.panels[1].discretized_elements[j]
                    coords_panel2 = elecoords
                    idx1p2 = idx1
                    idx2p2 = idx2
                    break
            if ele_panel1 is None or ele_panel2 is None:
                raise ValueError(
                    "Hinge line does not lie on the panels. Check node ordering.")

            possible_panel_ids1 = list(range(len(coords_panel1)))
            possible_panel_ids2 = list(range(len(coords_panel2)))
            possible_panel_ids1.remove(idx1p1)
            possible_panel_ids1.remove(idx2p1)
            possible_panel_ids2.remove(idx1p2)
            possible_panel_ids2.remove(idx2p2)
            nodepanel1 = possible_panel_ids1.pop()
            nodepanel2 = possible_panel_ids2.pop()

            new_hinge = [ele_panel1[nodepanel1], ele_panel1[idx1p1],
                         ele_panel1[idx2p1], ele_panel2[nodepanel2]]
            self.discretized_elements.append(new_hinge)


class Panel(Element):
    def __init__(self, nodes, thickness=0.35):
        self.thickness = thickness
        super().__init__(nodes)

    def mesh(self, n) -> Tuple[np.ndarray, list]:
        pass


class RectangularPanel(Panel):
    def __init__(self, nodes, thickness=0.35):
        super().__init__(nodes, thickness)
        self.opensees_type = 'ASDShellQ4'

    def mesh(self, n) -> Tuple[np.ndarray, list]:
        v1, v2, v3, v4 = np.array(self.coords)
        nodes = []
        node_index = {}

        for i in range(n+1):
            for j in range(n+1):
                point = (i/n)*v1 + (j/n)*v2 + ((n-i)/n)*v4 + ((n-j)/n)*v3
                idx = len(nodes)
                nodes.append(point)
                node_index[(i, j)] = idx

        elements = []
        for i in range(n):
            for j in range(n):
                a = (i, j)
                b = (i+1, j)
                c = (i+1, j+1)
                d = (i, j+1)
                elements.append(
                    [node_index[a], node_index[b], node_index[c], node_index[d]])
        new_nodes = np.array(nodes)
        new_elements = elements
        self.discretized_nodes = new_nodes
        self.discretized_elements = new_elements
        self._discretized_elements = [x[:] for x in new_elements]
        return new_nodes, new_elements


class TriangularPanel(Panel):
    def __init__(self, nodes, thickness=0.35):
        super().__init__(nodes, thickness)
        self.opensees_type = 'ASDShellT3'

    def mesh(self, n) -> Tuple[np.ndarray, list]:
        v1, v2, v3 = np.array(self.coords)
        nodes = []
        node_index = {}

        for i in range(n+1):
            for j in range(n+1-i):
                k = n - i - j
                point = (i/n)*v1 + (j/n)*v2 + (k/n)*v3
                idx = len(nodes)
                nodes.append(point)
                node_index[(i, j, k)] = idx

        elements = []
        for i in range(n):
            for j in range(n-i):
                k = n - i - j
                a = (i, j, k)
                b = (i+1, j, k-1)
                c = (i, j+1, k-1)
                d = (i+1, j+1, k-2)
                if k > 0:
                    elements.append(
                        [node_index[a], node_index[b], node_index[c]])
                if k > 1:
                    elements.append(
                        [node_index[b], node_index[d], node_index[c]])
        new_nodes = np.array(nodes)
        new_elements = elements
        self.discretized_nodes = new_nodes
        self.discretized_elements = new_elements
        self._discretized_elements = [x[:] for x in new_elements]
        return new_nodes, new_elements
