from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt


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
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
model.visualize(ax=ax, undeformed=True, node_labels=True)
plt.show()
a = 0

# self.geometry.tie_nodes

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
