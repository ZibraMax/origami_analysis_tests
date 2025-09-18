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
            panel.eletags = []
            for shell in panel.discretized_elements:
                nel = len(ops.getEleTags())
                panel.eletags.append(nel)
                ops.element(panel_type, nel, *shell,
                            self.materials_panels[i], '-corotational', '-reducedIntegration', '-drillingNL')

        for i, hinge in enumerate(self.geometry.hinges):
            if self.materials_hinges[i] == -1:
                raise ValueError(f"No material assigned to hinge {i}.")
            hinge.eletags = []
            for minihinge in hinge.discretized_elements:
                nel = len(ops.getEleTags())
                ops.element('OriHinge', nel, *minihinge,
                            self.materials_hinges[i], hinge.theta1, hinge.theta2)
                hinge.eletags.append(nel)
        for i, bar in enumerate(self.geometry.bars):
            nel = len(ops.getEleTags())
            bar.eletag = nel
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

        for stie in self.geometry.super_tie_nodes:
            ops.equalDOF(stie[0], stie[1], *[1, 2, 3, 4, 5, 6])

    def setup_model(self, tol=1e-5):
        ops.system('BandGeneral')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.test('NormDispIncr', tol, 500)

    def analyze(self, n_steps, callback=None):
        solutions = []
        try:
            for i in range(n_steps):
                if i != 0:
                    if ops.analyze(1) != 0:
                        print(f"Analysis failed at step {i}")
                        break
                if callback is not None:
                    res = callback(i)
                solutions.append(self.get_disp_vector(res))
        except KeyboardInterrupt as e:
            print(f"Analysis stopped at step {i} with error: {e}")
        self.solutions = solutions
        return solutions

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

        for bar in self.geometry.bars:
            n1, n2 = bar.discretized_nodes
            c1 = np.array(ops.nodeCoord(n1))
            c2 = np.array(ops.nodeCoord(n2))
            color = 'blue'
            alpha = 1
            disp1 = ops.nodeDisp(n1)
            disp2 = ops.nodeDisp(n2)
            c1 += np.array(disp1)[:3]*(1-undeformed)
            c2 += np.array(disp2)[:3]*(1-undeformed)
            ax.plot([c1[0], c2[0]], [c1[1], c2[1]], [c1[2], c2[2]],
                    color=color, alpha=alpha, lw=3)

        if plot_hinges:
            for hinge in self.geometry.hinges:
                for ele in hinge.discretized_elements:
                    n1, n2, n3, n4 = ele
                    c1 = np.array(ops.nodeCoord(n1)) + \
                        np.array(ops.nodeDisp(n1))[:3]*(1-undeformed)
                    c2 = np.array(ops.nodeCoord(n2)) + \
                        np.array(ops.nodeDisp(n2))[:3]*(1-undeformed)
                    c3 = np.array(ops.nodeCoord(n3)) + \
                        np.array(ops.nodeDisp(n3))[:3]*(1-undeformed)
                    c4 = np.array(ops.nodeCoord(n4)) + \
                        np.array(ops.nodeDisp(n4))[:3]*(1-undeformed)
                    color = 'red'
                    alpha = 1
                    # plot dash line between n1 and n2
                    ax.plot([c1[0], c2[0]], [c1[1], c2[1]], [c1[2], c2[2]],
                            color=color, alpha=alpha, linestyle='dashed', lw=3)
                    # plot solid line between n2 and n3
                    ax.plot([c2[0], c3[0]], [c2[1], c3[1]], [c2[2], c3[2]],
                            color=color, alpha=alpha, lw=3)
                    # plot dash line between n3 and n4
                    ax.plot([c3[0], c4[0]], [c3[1], c4[1]], [c3[2], c4[2]],
                            color=color, alpha=alpha, linestyle='dashed', lw=3)

        for i in ops.getNodeTags():
            coor = ops.nodeCoord(i) + np.array(ops.nodeDisp(i))[:3]
            color = 'yellow'
            size = 2
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

    def get_disp_vector(self, additional_info=None):
        additional_info = additional_info or {}
        U_NODES = []
        nodes = ops.getNodeTags()
        lam = ops.getLoadFactor(1)
        for j in nodes:
            U_NODES.extend(ops.nodeDisp(j))
        U_NODES = np.array(U_NODES)
        U_NODES = U_NODES.reshape((self.geometry.ngdl_per_node*len(nodes), 1))
        # Get integrator name
        return {
            "info": {**additional_info, "ld": lam},
            "U": U_NODES.tolist(),
        }

    def export_json(self, filename):
        data = self.geometry.to_json()

        materials_panels = []
        for i, mat in enumerate(self.materials_panels):
            materials_panels.extend(
                [mat]*len(self.geometry.panels[i].discretized_elements))

        data["properties"]["materials_panels"] = materials_panels
        data["properties"]["materials_bars"] = self.materials_bars
        data["properties"]["materials_hinges"] = self.materials_hinges
        data["regions"] = []
        data["ebc"] = []
        data["nbc"] = []
        data["nvn"] = 6
        if hasattr(self, "solutions"):
            data["solutions"] = self.solutions
        import json
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
