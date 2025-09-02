"""
Este ejemplo funciona con barras y hinges. Pero es demasiado sensible a los parametros del solver
"""


import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib.animation as animation


def visualize(ax=None, undeformed=False, node_labels=False):
    if ax is None:
        fig = plt.figure(figsize=[12, 5])
        ax = fig.add_subplot(1, 2, 1, projection='3d')
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
        if etype == "OriHinge" or etype == "BendHinge":
            ax.plot(coords[1:-1, 0], coords[1:-1, 1], coords[1:-1, 2], '--',
                    color='yellow', linewidth=2, zorder=100)
        elif etype == "CorotTruss":
            color = 'blue'
            if undeformed:
                color = 'black'
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
                    color=color)
        else:
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


data = json.load(open('./bar_and_hinge/instability/instability.json'))

nodes = data['nodes']
elements = data['dictionary']
types = data['types']
bars = [elements[i] for i in range(len(elements)) if types[i] == "CorotTruss"]
folds = [elements[i] for i in range(len(elements)) if types[i] == "OriHinge"]
bends = [elements[i] for i in range(len(elements)) if types[i] == "BendHinge"]
folds = folds[::2]

for i, element in enumerate(bends):
    bars.append([element[1], element[2]])
    types.append("CorotTruss")


types = ["CorotTruss"] * len(bars) + ["OriHinge"] * \
    len(folds) + ["BendHinge"] * len(bends)
elements = bars + folds + bends


ebc = data['ebc']
nbc = data['nbc']

nodes = np.array(nodes)
tolerance = 1e-6
fixed_nodes = np.where((np.abs(nodes[:, 0]) < tolerance) & (
    np.abs(nodes[:, 2]) < tolerance))[0]
for n in fixed_nodes:
    if not any(bv[0] == n for bv in ebc):
        ebc.append([int(n), 1, 1, 1])
# Now at x = maxx
maxx = np.max(nodes[:, 0])
fixed_nodes = np.where((np.abs(nodes[:, 0] - maxx) < tolerance) & (
    np.abs(nodes[:, 2]) < tolerance))[0]
for n in fixed_nodes:
    if not any(bv[0] == n for bv in ebc):
        ebc.append([int(n), 1, 1, 1])

minx = np.min(nodes[:, 0])
maxy = np.max(nodes[:, 1])
miny = np.min(nodes[:, 1])
center = [(maxx + minx) / 2, (maxy + miny) / 2-7.559289460184545/2]
distances = np.linalg.norm(nodes[:, 0:2] - center, axis=1)
n_load = np.argmin(distances)
nbc.append([int(n_load), 0, 0, -1.0])

E = 1000
A = 1
stifness_f = 10
k_hinges = 1
k_panels = stifness_f*k_hinges
fact = 1


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)
ops.uniaxialMaterial('Elastic', 1, E)

for i, v in enumerate(nodes):
    ops.node(i, *v)

for i, ele in enumerate(elements):
    nel = len(ops.getEleTags())
    if types[i] == "CorotTruss":
        ops.element('corotTruss', nel, *ele, A, 1)
    elif types[i] == "OriHinge":
        ops.element('OriHinge', nel, *ele, k_hinges, R(10), R(350))
    elif types[i] == "BendHinge":
        ops.element('OriHinge', nel, *ele, k_panels, R(10), R(350))

bcarray = np.array(ebc)
nodes_bc = bcarray[:, 0]
for bc in ebc:
    ops.fix(*bc)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

for load in nbc:
    ops.load(*load)

ops.system('Umfpack')
ops.constraints('Plain')
ops.numberer('RCM')
ops.algorithm('Newton')
ops.test('NormDispIncr', 1.0e-3, 100)
# ops.integrator('LoadControl', 0.1)
ops.integrator('MGDCM', 0.0005, 4, 2, 1)
ops.analysis('Static')
fig = plt.figure(figsize=[12, 5])
ax = fig.add_subplot(1, 1, 1, projection='3d')
visualize(ax=ax, undeformed=True, node_labels=False)
plt.show()

data = {}
data['step'] = []
data['load_factor'] = []
data['disp'] = []
M = 5000
n = 3


plt.ion()
fig, ax = plt.subplots()
line, = ax.plot([], [], 'r-')  # Initialize an empty plot
ax.grid()
ax.set_xlabel('Displacement')
ax.set_ylabel('Load factor')


def draw():
    line.set_data(data['disp'], data['load_factor'])
    ax.set_xlim(min(data['disp']), max(data['disp']) + 1)
    ax.set_ylim(min(data['load_factor']) - 1, max(data['load_factor']) + 1)
    plt.draw()
    plt.pause(0.1)  # Pause for a short duration to allow plot to update


flag = True
try:
    for i in range(1, M+1):
        if ops.analyze(1) != 0:
            print(f"Analysis failed at step {i}")
            break
        lam = ops.getLoadFactor(1)
        # theta = ops.eleResponse(len(dictionary)-2, 'theta')
        disp_z = ops.nodeDisp(nbc[0][0], 3)
        if i % 10 == 0:
            draw()

        print(f"{i},{lam},{disp_z}")
        data['step'].append(i)
        data['load_factor'].append(lam)
        data['disp'].append(-disp_z)
except KeyboardInterrupt as e:
    print(f"Analysis stopped at step {i} with error: {e}")

plt.ioff()
plt.show()
fig = plt.figure(figsize=[12, 5])
ax = fig.add_subplot(1, 2, 1, projection='3d')
visualize(plt.gca())
ax = plt.gcf().add_subplot(1, 2, 2)
plt.plot(data['disp'], data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Disp')
plt.ylabel('Load Factor')
plt.tight_layout()
plt.show()
