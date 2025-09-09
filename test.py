from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt


def visualize(ax=None, undeformed=False, node_labels=False, plot_hinges=True):
    if ax is None:
        fig = plt.figure(figsize=[12, 5])
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        plt.axis("off")

    etags = ops.getEleTags()
    dictionary = [ops.eleNodes(i) for i in etags]
    for i in etags:
        etype = ELE_TYPES[i]
        coords = np.array([ops.nodeCoord(n) for n in dictionary[i]])
        for j, n in enumerate(dictionary[i]):
            disp = ops.nodeDisp(n)
            coords[j, 0] += disp[0]
            coords[j, 1] += disp[1]
            coords[j, 2] += disp[2]
        if etype == "OriHinge" or etype == "BendHinge":
            if plot_hinges:
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
    for i in ops.getNodeTags():
        coor = ops.nodeCoord(i) + np.array(ops.nodeDisp(i))[:3]
        color = 'yellow'
        size = 0
        ax.scatter(coor[0], coor[1], coor[2], color=color, s=size)
        if node_labels:
            ax.text(coor[0], coor[1], coor[2], str(i), color='red', fontsize=8)
    if undeformed:
        # # Plot loads as unit lengt arrows in plot coords
        # for node, load in LOADS.items():
        #     coor = np.array(ops.nodeCoord(node))
        #     load = np.array(load)
        #     if np.linalg.norm(load) > 0:
        #         load = load / np.linalg.norm(load) * 0.1 * np.linalg.norm(np.array(ops.nodeCoord(
        #             ops.getNodeTags()[-1])) - np.array(ops.nodeCoord(ops.getNodeTags()[0])))
        #         coor = coor - load
        #         ax.quiver(coor[0], coor[1], coor[2], load[0], load[1],
        #                   load[2], color='red', arrow_length_ratio=0.5, linewidth=2)

        ax.axis('off')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')


H = 70
n = 6
b = 68

data = Kresling(b=b, H0=0.0, H1=H, n=n).generate()
O = Geometry.from_json(data, t=0.2)
O.mesh(n=2)

O.add_bc_plane('z', 0.0, values=[1, 1, 1, 0, 0, 0])
O.add_load_plane('z', H, values=[0, 0, -1/6, 0, 0, 0])

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3)

model = ShellAndHinge(O)
model.add_material_shells(mat_tag=1, E=100, v=0.3)
model.add_material_bars(mat_tag=2, E=1e9, A=1.0)
model.add_material_hinges(k=0.00)
model.create_model()
for n in nodes_tie[1:]:
    ops.equalDOF(nodes_tie[0], n, 3)


ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Plain')
ops.test('NormDispIncr', 1.0e-3, 500)

M = 500
ops.integrator('DisplacementControl', nodes_tie[0], 3, -70/M)
ops.algorithm('Newton')
ops.analysis('Static')


fig = plt.figure(figsize=[12, 6])
ax = fig.add_subplot(1, 1, 1)
N = 2
res = {"step": [], "load_factor": [], "disp": []}
flag = True
try:
    for i in range(0, M+1):
        if i != 0:
            if ops.analyze(1) != 0:
                print(f"Analysis failed at step {i}")
                break
        lam = ops.getLoadFactor(1)
        disp_z = ops.nodeDisp(nodes_tie[0], 3)
        print(f"{i},{lam},{disp_z}")
        res['step'].append(i)
        res['load_factor'].append(lam)
        res['disp'].append(-disp_z)

except KeyboardInterrupt as e:
    print(f"Analysis stopped at step {i} with error: {e}")

ax.plot(res['disp'], res['load_factor'], 'r-')
ax.set_xlabel('Displacement')
ax.set_ylabel('Load factor')
ax.grid()

plt.tight_layout()
plt.show()
