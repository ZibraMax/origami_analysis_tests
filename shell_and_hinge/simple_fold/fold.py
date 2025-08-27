import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json


def visualize(ax=None, undeformed=False):
    if ax is None:
        fig = plt.figure(figsize=[12, 5])
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        # plt.axis('off')
        # Add arrow for load in 3D. Use a quiver plot with a single point
        ax.quiver(*[-1.1, 0.0, 0.4330127018922196], 0.23, 0, 0, color='k', linewidth=2,
                  arrow_length_ratio=0.5, zorder=300)
        ax.text(*[-1.15, 0.0, 0.4330127018922196],
                'P', color='k', fontsize=12)
        # Add small X,Y,Z axes
        ax.quiver(*[-0.5, -1.0, 0.0], 0.2, 0, 0, color='r', linewidth=1,
                  arrow_length_ratio=0.5, zorder=300)
        ax.quiver(*[-0.5, -1.0, 0.0], 0, 0.2, 0, color='g', linewidth=1,
                  arrow_length_ratio=0.5, zorder=300)
        ax.quiver(*[-0.5, -1.0, 0.0], 0, 0, 0.2, color='b', linewidth=1,
                  arrow_length_ratio=0.5, zorder=300)
        ax.text(1.25-1.5, -1.0, 0.0, 'X', color='r', fontsize=12)
        ax.text(1.0-1.5, -0.75, 0.0, 'Y', color='g', fontsize=12)
        ax.text(1.0-1.5, -1.0, 0.25, 'Z', color='b', fontsize=12)
        plt.axis("off")

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
            ax.plot(coords[1:-1, 0], coords[1:-1, 1], coords[1:-1, 2], '--',
                    color='yellow', linewidth=2, zorder=100)
        else:
            # Close the loop
            color = 'green'
            alpha = 0.5
            if undeformed:
                color = 'black'
                alpha = 1
            coords = np.vstack((coords, coords[0, :]))
            ax.plot_trisurf(coords[:, 0], coords[:, 1], coords[:, 2],
                            color=color, alpha=alpha)
    # Plot nodes as points
    for i in range(len(nodes)):
        coor = ops.nodeCoord(i)
        ax.scatter(coor[0], coor[1], coor[2], color='yellow', s=20, zorder=200)

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

EXTRA_NODES = True
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
        ops.element('ASDShellT3', i, *element, 1, '-corotational')

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
ops.test('NormDispIncr', 1.0e-3, 200)
ops.algorithm('Newton')
ops.integrator('ArcLength', 0.1, 0.1)
# ops.integrator('LoadControl', 0.1)
ops.analysis('Static')


data = {}
data['step'] = []
data['load_factor'] = []
data['angle_degrees'] = []
data['disp_node4_z'] = []

print("Step,Load_Factor,Angle(degrees),Disp_Node4_Z")
M = 60
n = 10
visualize(undeformed=True)
for i in range(1, M+1):
    if ops.analyze(1) != 0:
        print(f"Analysis failed at step {i}")
        break
    lam = ops.getLoadFactor(1)
    theta = ops.eleResponse(2, 'theta')
    disp_z = ops.nodeDisp(3, 3)
    print(f"{i},{lam},{theta[0] * 180.0 / np.pi},{disp_z}")
    data['step'].append(i)
    data['load_factor'].append(lam)
    data['angle_degrees'].append(theta[0] * 180.0 / np.pi)
    data['disp_node4_z'].append(disp_z)
    if i % (M//n) == 0 or i == M:
        visualize(plt.gca())
ax = plt.gcf().add_subplot(1, 2, 2)
plt.plot(data['angle_degrees'], data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Angle (degrees)')
plt.ylabel('Load Factor')
plt.tight_layout()
plt.show()
