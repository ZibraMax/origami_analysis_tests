from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt


H = 70
H0 = 35
n = 6
number_sides = n
b = 52
khinge = 0.02
thickness = 0.5
kresling = Kresling(b=b, H0=H0, H1=H, n=n)
data = kresling.generate(get_int_lines=False,
                         get_ext_lines=True,
                         get_base_bars=True,
                         get_ext_hinges=False,
                         get_int_hinges=True,
                         get_panels=True,
                         get_base_panels=True,
                         get_base_hinges=True)
delta_theta = kresling.delta_theta*1.08
data['properties']['theta1'] = R(20)
data['properties']['theta2'] = R(340)
O = Geometry.from_json(data, t=thickness)
O.mesh(n=5, mesh_hinges=True)
hinge = O.hinges[-1]
res, center_idx = O.test_point_vertices([0.0, 0.0, H], tol=1e-3)

O.add_bc_plane('z', 0.0, values=[1, 1, 1, 0, 0, 0])

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3, basenodes=True)
print("Nodes at the top plane:", nodes_tie)

# O.add_load_plane('z', H, values=[0, 0, -1/len(nodes_tie), 0, 0, 0])
O.nbc.append([center_idx, 0, 0, -1, 0, 0, 0])

BASE_E = 210000
model = ShellAndHinge(O)
model.add_material_shells(mat_tag=1, E=BASE_E, v=0.49)
model.add_material_shells(mat_tag=2, E=1e6*BASE_E, v=0.49,
                          shell_list=list(range(2*n, 2*n+2)))
model.add_material_bars(mat_tag=3, E=1e10, A=1.0)
model.add_material_hinges(k=khinge)
model.create_model()

nodes_tie = O.get_nodes_plane('z', H, tol=1e-3)
for n in nodes_tie[1:]:
    ops.equalDOF(nodes_tie[0], n, 3)

ops.fix(*[center_idx, 1, 1, 0, 0, 0, 0])


# get centernode


model.setup_model(tol=1e-3)

# Setup solver
M = 100
ops.integrator('DisplacementControl', center_idx, 3, -H/M)
# ops.integrator('MGDCM', 200, 15, 4, 0)
ops.algorithm('Newton')
ops.analysis('Static')

Nmodes = 12
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
    disp_z = ops.nodeDisp(center_idx, 3)
    rot_top = ops.nodeDisp(center_idx, 6)
    print(f"{i},{lam},{disp_z}")
    res['step'].append(i)
    res['load_factor'].append(lam)
    res['disp'].append(-disp_z)

    print(",".join([f"{D(i)}" for i in hinge.get_theta()]))
    return {"solver-type": "DisplacementControl", "vertical-disp": disp_z, "rotation_top_node": rot_top}


solutions = model.analyze(M, callback=callback)

fig = plt.figure(figsize=[12, 6])
ax2 = fig.add_subplot(1, 2, 1, projection='3d')
ax = fig.add_subplot(1, 2, 2)
model.visualize(ax=ax2)
model.export_json(
    f"kresling_n_{number_sides}_b_{b}_h_{H}_h0_{H0}_t_{thickness}_kf_{khinge}.json")
ax.plot(res['disp'], res['load_factor'], 'r-')
ax.set_xlabel('Displacement')
ax.set_ylabel('Load factor')
ax.grid()
plt.tight_layout()
plt.show()
