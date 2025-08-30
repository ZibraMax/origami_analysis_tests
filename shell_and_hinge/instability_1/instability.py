import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json


def create_extra_nodes_for_hinges(nodes, elemements, types):
    if not isinstance(nodes, np.ndarray):
        nodes = np.array(nodes)
    nodes = nodes.tolist()
    nodes_hinges_middle = []
    nodes_to_element = {}
    tie_nodes = []
    for i, etype in enumerate(types):
        if etype == "OriHinge":
            n1, n2, n3, n4 = elemements[i]
            nodes_hinges_middle.append(n2)
            nodes_hinges_middle.append(n3)
        else:
            element = elemements[i]
            for n in element:
                if n not in nodes_to_element:
                    nodes_to_element[n] = []
                nodes_to_element[n].append(i)
    nodes_hinges_middle = list(set(nodes_hinges_middle))
    for node in nodes_hinges_middle:
        connected_elements = nodes_to_element[node]
        number_conected = len(connected_elements)
        for i in range(number_conected - 1):
            new_node_id = len(nodes)
            nodes.append(nodes[node])
            tie_nodes.append((node, new_node_id))
            elem_idx = connected_elements[i+1]
            pos_new_node = elemements[elem_idx].index(node)
            elemements[elem_idx][pos_new_node] = new_node_id

    return np.array(nodes), elemements, tie_nodes


def visualize(ax=None, undeformed=False, node_labels=False):
    if ax is None:
        fig = plt.figure(figsize=[12, 5])
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        # plt.axis('off')
        # Add arrow for load in 3D. Use a quiver plot with a single point
        # ax.quiver(*[-1.1, 0.0, 0.4330127018922196], 0.23, 0, 0, color='k', linewidth=2,
        #           arrow_length_ratio=0.5, zorder=300)
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
        coor = ops.nodeCoord(i) + np.array(ops.nodeDisp(i))[:3]
        # choose color based on whether node is constrained
        color = 'yellow'
        size = 0
        if any(bv[0] == i for bv in ebc):
            color = 'red'
            size = 20
        elif any(nc[0] == i for nc in nbc):
            color = 'blue'
            size = 20
        ax.scatter(coor[0], coor[1], coor[2], color=color, s=size)
        if node_labels:
            ax.text(coor[0], coor[1], coor[2], str(i), color='red', fontsize=8)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')


def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


THETA_1 = R(10)
THETA_2 = R(350)

data = json.load(open('./shell_and_hinge/instability_1/geometry.json'))
nodes = data['nodes']
dictionary = data['dictionary']
types = data['types']

ebc = data['ebc']
nbc = data['nbc']
props = data['properties']
E = props['E']
v = props['v']
t = props['t']
kf = props['kf']

# detect nodes at x = 0 and z =0
nodes = np.array(nodes)
tolerance = 1e-6
fixed_nodes = np.where((np.abs(nodes[:, 0]) < tolerance) & (
    np.abs(nodes[:, 2]) < tolerance))[0]
for n in fixed_nodes:
    if not any(bv[0] == n for bv in ebc):
        ebc.append([int(n), 1, 1, 1, 0, 0, 0])
# Now at x = maxx
maxx = np.max(nodes[:, 0])
fixed_nodes = np.where((np.abs(nodes[:, 0] - maxx) < tolerance) & (
    np.abs(nodes[:, 2]) < tolerance))[0]
for n in fixed_nodes:
    if not any(bv[0] == n for bv in ebc):
        ebc.append([int(n), 1, 1, 1, 0, 0, 0])


# Add a single load to the closest node to the center of the domain
minx = np.min(nodes[:, 0])
maxy = np.max(nodes[:, 1])
miny = np.min(nodes[:, 1])
center = [(maxx + minx) / 2, (maxy + miny) / 2-7.559289460184545/2]
distances = np.linalg.norm(nodes[:, 0:2] - center, axis=1)
n_load = np.argmin(distances)
nbc.append([int(n_load), 0, 0, -1.0, 0, 0, 0])

dictionary_shells = [i for i, etype in enumerate(
    data['types']) if etype != "OriHinge"]
dictionary_hinges = [i for i, etype in enumerate(
    data['types']) if etype == "OriHinge"]
dictionary = np.array(dictionary).astype(
    int)[dictionary_shells + dictionary_hinges[::2]].tolist()
types = np.array(types)[dictionary_shells +
                        dictionary_hinges[::2]].tolist()

EXTRA_NODES = True
if EXTRA_NODES:
    nodes, dictionary, pairs_constraints = create_extra_nodes_for_hinges(
        nodes, dictionary, types)

nodes = np.array(nodes)


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
    elif etype == 'ASDShellQ4':
        ops.element('ASDShellQ4', i, *element, 1, '-corotational')

# Apply boundary conditions
for bv in ebc:
    ops.fix(*bv)
# Apply loads
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
for load in nbc:
    ops.load(*load)


ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('NormDispIncr', 1.0e-3, 100)
# ops.integrator('LoadControl', 0.1)
ops.integrator('DisplacementControl', nbc[0][0], 3, -0.1)
ops.algorithm('Newton')
ops.analysis('Static')

data = {}
data['step'] = []
data['load_factor'] = []
data['disp'] = []

print("Step,Load_Factor")
M = 5000
n = 3
fig = plt.figure(figsize=[12, 5])
ax = fig.add_subplot(projection='3d')
visualize(ax=ax, undeformed=False, node_labels=True)
plt.tight_layout()
plt.show()

fig = plt.figure(figsize=[12, 5])
ax = fig.add_subplot(1, 2, 1, projection='3d')
try:
    for i in range(1, M+1):
        if ops.analyze(1) != 0:
            print(f"Analysis failed at step {i}")
            break
        lam = ops.getLoadFactor(1)
        # theta = ops.eleResponse(len(dictionary)-2, 'theta')
        disp_z = ops.nodeDisp(nbc[0][0], 3)
        if abs(disp_z) > 35:
            print(f"Switching to arclength {i}")
            ops.integrator('ArcLength', 0.01, 0.01)

        print(f"{i},{lam},{disp_z}")
        data['step'].append(i)
        data['load_factor'].append(lam)
        data['disp'].append(-disp_z)
except KeyboardInterrupt as e:
    print(f"Analysis stopped at step {i} with error: {e}")
visualize(plt.gca())
ax = plt.gcf().add_subplot(1, 2, 2)
plt.plot(data['disp'], data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Disp')
plt.ylabel('Load Factor')
plt.tight_layout()
plt.show()
