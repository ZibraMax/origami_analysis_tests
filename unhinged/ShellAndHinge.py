import numpy as np
import opensees as ops
from .Geometry import Geometry, RectangularPanel, TriangularPanel
import matplotlib.pyplot as plt


class ShellAndHinge():
    def __init__(self, geometry: Geometry):
        self.geometry = geometry
        self.ndm = 3
        self.ndf = geometry.ngdl_per_node
        self.materials_panels = []
        self.materials_bars = []
        self.materials_hinges = []

        self.uniques_t = []
        for panel in self.geometry.panels:
            if panel.thickness not in self.uniques_t:
                self.uniques_t.append(panel.thickness)
        if len(self.uniques_t) > 1:
            print(
                "Warning: Multiple thicknesses detected. This is not supported yet.")
        ops.wipe()
        ops.model('basic', '-ndm', self.ndm, '-ndf', self.ndf)

    def add_material_shells(self, mat_tag, E, v, shell_list=None):
        if shell_list is None:
            self.materials_panels = [mat_tag]*len(self.geometry.panels)
        else:
            for i in shell_list:
                if i < 0 or i >= len(self.geometry.panels):
                    raise ValueError(f"Panel index {i} out of range.")
            for i in shell_list:
                while len(self.materials_panels) <= i:
                    self.materials_panels.append(-1)
                self.materials_panels[i] = mat_tag
        ops.section('ElasticMembranePlateSection',
                    mat_tag, E, v, self.uniques_t[0])

    def add_material_bars(self, mat_tag, E, A, bar_list=None):
        if bar_list is None:
            self.materials_bars = [mat_tag]*len(self.geometry.bars)
        else:
            for i in bar_list:
                if i < 0 or i >= len(self.geometry.bars):
                    raise ValueError(f"Bar index {i} out of range.")
            for i in bar_list:
                while len(self.materials_bars) <= i:
                    self.materials_bars.append(-1)
                self.materials_bars[i] = mat_tag
        ops.uniaxialMaterial('Elastic', mat_tag, E*A)

    def add_material_hinges(self, k, hinge_list=None):
        if hinge_list is None:
            self.materials_hinges = [k]*len(self.geometry.hinges)
        else:
            for i in hinge_list:
                if i < 0 or i >= len(self.geometry.hinges):
                    raise ValueError(f"Hinge index {i} out of range.")
            for i in hinge_list:
                while len(self.materials_hinges) <= i:
                    self.materials_hinges.append(-1)
                self.materials_hinges[i] = k

    def create_model(self):
        self.create_nodes()
        self.create_elements()
        self.set_ebcs()
        self.set_nbc()
        self.create_ties()

    def create_nodes(self):
        for i, node in enumerate(self.geometry.nodes):
            ops.node(i, *node)

    def create_elements(self):
        for i, panel in enumerate(self.geometry.panels):
            if self.materials_panels[i] == -1:
                raise ValueError(f"No material assigned to panel {i}.")

            panel_type = panel.opensees_type
            for shell in panel.discretized_elements:
                nel = len(ops.getEleTags())
                ops.element(panel_type, nel, *shell,
                            self.materials_panels[i], '-corotational')

        for i, hinge in enumerate(self.geometry.hinges):
            nel = len(ops.getEleTags())
            if self.materials_hinges[i] == -1:
                raise ValueError(f"No material assigned to hinge {i}.")
            ops.element('OriHinge', nel, *hinge.discretized_nodes,
                        self.materials_hinges[i], hinge.theta1, hinge.theta2)
        for i, bar in enumerate(self.geometry.bars):
            nel = len(ops.getEleTags())
            if self.materials_bars[i] == -1:
                raise ValueError(f"No material assigned to bar {i}.")
            ops.element('CorotTruss', nel, *bar.discretized_nodes, 1,
                        self.materials_bars[i])

    def set_ebcs(self):
        for ebc in self.geometry.ebc:
            ops.fix(*ebc)

    def set_nbc(self):
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        for nbc in self.geometry.nbc:
            ops.load(*nbc)

    def create_ties(self, tie_type=[1, 2, 3]):
        for tie in self.geometry.tie_nodes:
            ops.equalDOF(tie[0], tie[1], *tie_type)

    def visualize(self, ax=None, undeformed=False, node_labels=False, plot_hinges=True):
        if ax is None:
            fig = plt.figure(figsize=[12, 5])
            ax = fig.add_subplot(1, 2, 1, projection='3d')
            plt.axis("off")

        for panel in self.geometry.panels:
            for shell in panel.discretized_elements:

                coords = np.array([self.geometry.nodes[n] for n in shell])
                coords = coords.copy()
                for j, n in enumerate(shell):
                    disp = ops.nodeDisp(n)
                    coords[j, 0] += disp[0]*(1-undeformed)
                    coords[j, 1] += disp[1]*(1-undeformed)
                    coords[j, 2] += disp[2]*(1-undeformed)
                color = 'green'
                alpha = 0.5
                if undeformed:
                    color = 'black'
                    alpha = 0.8
                coords = np.vstack((coords, coords[0, :]))
                ax.plot_trisurf(coords[:, 0], coords[:, 1], coords[:, 2],
                                color=color, alpha=alpha)
        for i in ops.getNodeTags():
            coor = ops.nodeCoord(i) + np.array(ops.nodeDisp(i))[:3]
            color = 'yellow'
            size = 0
            ax.scatter(coor[0], coor[1], coor[2], color=color, s=size)
            if node_labels:
                ax.text(coor[0], coor[1], coor[2],
                        str(i), color='red', fontsize=8)
        for tie in self.geometry.tie_nodes:
            c1 = ops.nodeCoord(tie[0])
            c2 = ops.nodeCoord(tie[0])
            if c1 != c2:
                print("Tie nodes are not coincident!")
            ax.scatter(*c1, color='yellow', s=10)
        if undeformed:
            # Plot loads as unit lengt arrows in plot coords
            for a in self.geometry.nbc:
                node, load = a[0], a[1:]
                load = np.array(load)[:3]
                coor = np.array(ops.nodeCoord(node))
                load = np.array(load)
                if np.linalg.norm(load) > 0:
                    load = load / np.linalg.norm(load) * 0.1 * np.linalg.norm(np.array(ops.nodeCoord(
                        ops.getNodeTags()[-1])) - np.array(ops.nodeCoord(ops.getNodeTags()[0])))
                    coor = coor - load
                    ax.quiver(coor[0], coor[1], coor[2], load[0], load[1],
                              load[2], color='red', arrow_length_ratio=0.5, linewidth=2)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_aspect('equal')

    def get_disp_vector(self):
        U_NODES = []
        nodes = ops.getNodeTags()
        lam = ops.getLoadFactor(1)
        for j in nodes:
            U_NODES.extend(ops.nodeDisp(j))
        U_NODES = np.array(U_NODES)
        U_NODES = U_NODES.reshape((self.geometry.ngdl_per_node*len(nodes), 1))
        # Get integrator name
        return {
            "info": {"solver-type": "int_name", "ld": lam},
            "U": U_NODES.tolist(),
        }
