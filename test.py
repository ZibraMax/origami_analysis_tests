from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt


H = 100
n = 6
b = 68

data = Kresling(b=b, H0=0, H1=H, n=n).generate(get_int_lines=False)

O = Geometry.from_json(data, t=0.2)
O.mesh(n=2)

O.add_bc_plane('z', 0.0, values=[1, 1, 1, 0, 0, 0])

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3)
print("Nodes at the top plane:", nodes_tie)

O.add_load_plane('z', H, values=[0, 0, -1/len(nodes_tie), 0, 0, 0])

model = ShellAndHinge(O)
model.add_material_shells(mat_tag=1, E=100, v=0.3)
model.add_material_bars(mat_tag=2, E=1e9, A=1.0)
model.add_material_hinges(k=0.00)
model.create_model()

# Extra manual ties
nodes_tie = O.get_nodes_plane('z', H, tol=1e-3)
for n in nodes_tie[1:]:
    ops.equalDOF(O.node_map[nodes_tie[0]][0], O.node_map[n][0], 3)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
model.visualize(ax=ax, undeformed=True, node_labels=True)
plt.show()

model.setup_model()

# Setup solver
M = 500
ops.integrator('DisplacementControl', nodes_tie[0], 3, -H/M)
ops.algorithm('Newton')
ops.analysis('Static')


res = {"step": [], "load_factor": [], "disp": []}


def callback(i):
    lam = ops.getLoadFactor(1)
    disp_z = ops.nodeDisp(nodes_tie[0], 3)
    print(f"{i},{lam},{disp_z}")
    res['step'].append(i)
    res['load_factor'].append(lam)
    res['disp'].append(-disp_z)


fig = plt.figure(figsize=[12, 6])
ax2 = fig.add_subplot(1, 2, 1, projection='3d')
ax = fig.add_subplot(1, 2, 2)
solutions = model.analyze(M, callback=callback)
model.visualize(ax=ax2)
model.export_json("out_kresling.json")
ax.plot(res['disp'], res['load_factor'], 'r-')
ax.set_xlabel('Displacement')
ax.set_ylabel('Load factor')
ax.grid()
plt.tight_layout()
plt.show()
