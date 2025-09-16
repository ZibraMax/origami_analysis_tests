from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt


H = 70
n = 6
b = 68
kresling = Kresling(b=b, H0=0.0, H1=H, n=n)
data = kresling.generate(get_int_lines=False,
                         get_ext_lines=True,
                         get_base_bars=False,
                         get_ext_hinges=False,
                         get_int_hinges=True,
                         get_panels=True,
                         get_base_panels=True)
delta_theta = kresling.delta_theta
O = Geometry.from_json(data, t=1.5)
O.mesh(n=4)

O.add_bc_plane('z', 0.0, values=[1, 1, 1, 0, 0, 0])

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3, basenodes=True)
print("Nodes at the top plane:", nodes_tie)

O.add_load_plane('z', H, values=[0, 0, -1/len(nodes_tie), 0, 0, 0])

model = ShellAndHinge(O)
model.add_material_shells(mat_tag=1, E=10000000, v=0.3)
model.add_material_shells(mat_tag=2, E=40000000, v=0.3,
                          shell_list=list(range(2*n, 2*n+2)))
# model.add_material_bars(mat_tag=3, E=100, A=1.0)
model.add_material_hinges(k=0.2)
model.create_model()

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3)
for n in nodes_tie[1:]:
    ops.equalDOF(nodes_tie[0], n, 3)
res, center_idx = O.test_point_vertices([0.0, 0.0, H], tol=1e-3)

ops.fix(*[center_idx, 1, 1, 0, 0, 0, 0])


# get centernode


model.setup_model(tol=1e-3)

# Setup solver
M = 200
ops.integrator('DisplacementControl', center_idx, 6, -delta_theta/M)
# ops.integrator('MGDCM', 0.2, 15, 4, 0)
ops.algorithm('Newton')
ops.analysis('Static')

Nmodes = 4
lam = ops.eigen('standard', 'symmBandLapack', Nmodes)
eigenvectors = []
for node in ops.getNodeTags():
    eigenvectors.append([])
    for mode in range(Nmodes):
        ev = ops.nodeEigenvector(node, mode+1)
        eigenvectors[-1].append(ev)

print("First eigenvalues:", lam)
model.solutions = []
factors = [0 for i in lam]
for mode in range(Nmodes):
    for node in ops.getNodeTags():
        nodedisp = ops.nodeDisp(node)
        for i, d in enumerate(nodedisp):
            ops.setNodeDisp(
                node, i+1, factors[mode]*eigenvectors[node][mode][i], '-commit')
    sol = model.get_disp_vector()
    sol["info"] = {"solver-type": "EIGEN", "ld": lam[mode]}
    model.solutions.append(sol)

model.export_json("eigv_kresling2.json")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
model.visualize(ax=ax, undeformed=False, node_labels=False)
plt.show()


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
