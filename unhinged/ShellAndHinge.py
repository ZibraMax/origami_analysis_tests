import numpy as np
import opensees as ops
from .Geometry import Geometry, RectangularPanel, TriangularPanel


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
