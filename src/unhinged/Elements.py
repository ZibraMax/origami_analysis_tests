import numpy as np
import opensees as ops
from typing import Tuple
def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


class Element():
    def __init__(self, nodes):
        self.nodes = nodes
        self.coords = []
        self.gdl = []
        self.discretized_nodes = []
        self.eletag = None

    def set_domanin(self, domain):
        self.coords = domain.base_nodes[self.nodes]
        self.gdl = domain.base_gdls[self.nodes]


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

    def get_theta(self):
        return ops.elementResponse(self.eletag, 'theta')


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
        return new_nodes, new_elements
