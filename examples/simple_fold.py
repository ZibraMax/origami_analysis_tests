import unhinged as uh
from unhinged import R, D
import numpy as np
import matplotlib.pyplot as plt
import opensees as ops
nodes = [
    [0.0, -0.5, 0.0],
    [0.0, 0.5, 0.0],
    [0.8660254037844386, 0.0, 0.0],
    [-0.7499999999999998, 0.0, -0.4330127018922196]
]
kf = 0.3
E = 200
v = 0.3
t = 0.5
N = 3

geo = uh.Geometry(6)
geo.add_nodes(nodes)

panel1 = uh.TriangularPanel([0, 2, 1], thickness=t)
panel2 = uh.TriangularPanel([0, 1, 3], thickness=t)
hinge = uh.Hinge([3, 0, 1, 2])

geo.add_panel(panel1)
geo.add_panel(panel2)
geo.add_hinge(hinge)

geo.mesh(N)
geo.add_bc_plane('z', 0.0, [1, 1, 1, 0, 0, 0])
# geo.add_ebc(3, [1, 0, 0, 0, 0, 0])
geo.add_nbc(3, [1, 0, 0, 0, 0, 0])
model = uh.ShellAndHinge(geo)
model.add_material_hinges(kf)
model.add_material_shells(1, E, v)
model.create_model()
model.setup_model()


ops.integrator('MGDCM', 0.005, 15, 4, 0)
# ops.integrator('LoadControl', 0.01)
# ops.integrator('DisplacementControl', geo.node_map[3][0], 3, 0.005)
ops.algorithm('Newton')
ops.analysis('Static')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
model.visualize(ax=ax, undeformed=False, node_labels=False)
plt.show()

res = {"step": [], "load_factor": [], "disp": []}

M = 1000


def callback(i):
    lam = ops.getLoadFactor(1)
    disp_z = ops.nodeDisp(model.geometry.node_map[3][0], 3)
    print(f"{i},{lam},{disp_z}")
    res['step'].append(i)
    res['load_factor'].append(lam)
    res['disp'].append(disp_z)


fig = plt.figure(figsize=[12, 6])
ax2 = fig.add_subplot(1, 2, 1, projection='3d')
ax = fig.add_subplot(1, 2, 2)
solutions = model.analyze(M, callback=callback)
model.visualize(ax=ax2)
model.export_json("out_fold.json")
ax.plot(res['disp'], res['load_factor'], 'r-')
ax.set_xlabel('Displacement')
ax.set_ylabel('Load factor')
ax.grid()
plt.tight_layout()
plt.show()
