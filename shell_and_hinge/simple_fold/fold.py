import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json


def visualize(ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

    etags = ops.getEleTags()
    dictionary = [ops.eleNodes(i) for i in etags]
    for i in etags:
        etype = types[i]
        coords = np.array([ops.nodeCoord(n) for n in dictionary[i]])
        for j, n in enumerate(dictionary[i]):
            disp = ops.nodeDisp(n)
            coords[j, 0] += disp[0]
            coords[j, 1] += disp[1]
            coords[j, 2] += disp[2]
        if etype == "OriHinge":
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], '--',
                    color='yellow', linewidth=3, zorder=100)
        else:
            # Close the loop
            coords = np.vstack((coords, coords[0, :]))
            ax.plot_trisurf(coords[:, 0], coords[:, 1], coords[:, 2],
                            color='black', alpha=0.6)
    # Plot nodes as points
    for i in range(len(nodes)):
        coor = ops.nodeCoord(i)
        ax.scatter(coor[0], coor[1], coor[2], color='red', s=20, zorder=200)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')


def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


THETA_1 = R(10)
THETA_2 = R(350)

data = json.load(open('./shell_and_hinge/simple_fold/fold.json'))

nodes = data['nodes']
dictionary = data['dictionary']

dictionary_shells = [i for i, etype in enumerate(
    data['types']) if etype != "OriHinge"]
dictionary_hinges = [i for i, etype in enumerate(
    data['types']) if etype == "OriHinge"]
print(dictionary_hinges)
# For every hinge element, identify the two faces it connects
hinges_to_faces = {}
for i in dictionary_hinges:
    n1, n2, n3, n4 = data['dictionary'][i]
    face1 = None
    face2 = None
    for j in dictionary_shells:
        shell = data['dictionary'][j]
        if n1 in shell and n2 in shell and n3 in shell:
            face1 = j
        if n2 in shell and n3 in shell and n4 in shell:
            face2 = j
    hinges_to_faces[i] = (face1, face2)

pairs_constraints = []

EXTRA_NODES = False
if EXTRA_NODES:
    for i, faces in hinges_to_faces.items():
        face1, face2 = faces
        n1, n2, n3, n4 = data['dictionary'][i]
        shell = data['dictionary'][face2]

        idx_selected = shell.index(n2)
        nodes.append(nodes[n2])
        new_node = len(nodes)-1
        pairs_constraints.append((n2, new_node))
        dictionary[face2][idx_selected] = new_node

        idx_selected = shell.index(n3)
        nodes.append(nodes[n3])
        new_node = len(nodes)-1
        pairs_constraints.append((n3, new_node))
        dictionary[face2][idx_selected] = new_node

nodes = np.array(nodes)
types = data['types']
ebc = data['ebc']
nbc = data['nbc']
nvn = data['nvn']
props = data['properties']
E = props['E']
v = props['v']
t = props['t']
kf = props['kf']


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
ops.section('ElasticMembranePlateSection', 1, E, v, t)

# Add all nodes
for i, coord in enumerate(nodes):
    ops.node(i, *coord)

for pair in pairs_constraints:
    ops.equalDOF(pair[0], pair[1], 1, 2, 3)
# Iterate over elements (dictionary)
for i, element in enumerate(dictionary):
    etype = types[i]
    if etype == 'OriHinge':
        ops.element('OriHinge', i, *element, kf, THETA_1, THETA_2)
    elif etype == 'ASDShellT3':
        ops.element('ASDShellT3', i, *element, 1)

# Apply boundary conditions
for bv in ebc:
    ops.fix(*bv)
# Apply loads
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
for load in nbc:
    ops.load(*load)
ops.system('Umfpack')
ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-5, 100)
ops.algorithm('Newton')
# ops.integrator('ArcLength', 0.1, 0.1)
ops.integrator('LoadControl', 0.1)
ops.analysis('Static')

M = 100
n = 5
visualize()
for i in range(1, M+1):
    if ops.analyze(1) != 0:
        print(f"Analysis failed at step {i}")
        break
    lam = ops.getLoadFactor(1)
    theta = ops.eleResponse(2, 'theta') or [0]
    disp_z = ops.nodeDisp(3, 3)
    print(f"{i},{lam},{theta[0] * 180.0 / np.pi},{disp_z}")
    if i % (M//n) == 0 or i == M:
        visualize(plt.gca())
plt.show()
