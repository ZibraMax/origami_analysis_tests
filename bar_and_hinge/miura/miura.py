import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json


def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


data = json.load(open('./bar_and_hinge/miura/miura.json'))
# data has vertices, lines, perimeters, folds, bending, essential boundary condityions (ebc), loads (nbc)
vertices = data['vertices']
lines = data['lines']
folds = data['folds']
bending = data['bending']
ebc = data['ebc']
nbc = data['nbc']
E = 1000
A = 1
stifness_f = 100
k_hinges = 0.1
k_panels = stifness_f*k_hinges
fact = 1


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.uniaxialMaterial('Elastic', 1, E)


for i, v in enumerate(vertices):
    ops.node(i, *v)

for i, bars in enumerate(lines):
    nel = len(ops.getEleTags())
    ops.element('corotTruss', nel, *bars, A, 1)

hinges_element_numbers = []
for i, bending in enumerate(bending):
    nel = len(ops.getEleTags())
    n1, n2, n3, n4 = bending
    ops.element('OriHinge', nel, n1, n2, n3, n4, k_panels, R(10), R(350))
    hinges_element_numbers.append(nel)

for i, fold in enumerate(folds):
    nel = len(ops.getEleTags())
    n1, n2, n3, n4 = fold
    ops.element('OriHinge', nel, n1, n2, n3, n4, k_hinges, R(10), R(350))
    hinges_element_numbers.append(nel)

bcarray = np.array(ebc)
nodes_bc = bcarray[:, 0]

for node in range(len(vertices)):
    if node not in nodes_bc:
        ops.fix(node, 0, 0, 0, 1, 1, 1)
for bc in ebc:
    ops.fix(bc[0], *(bc[1:]+[1, 1, 1]))

ops.equalDOF(2, 5, 1)
ops.equalDOF(5, 8, 1)


ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
for load in nbc:
    ops.load(load[0], *np.array(load[1:]+[0, 0, 0])*fact)
ops.system('Umfpack')
ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormUnbalance', 1e-5, 50)
ops.algorithm('Newton')
# ops.integrator('ArcLength', 1, 1)
# ops.integrator('LoadControl', 0.01)
ops.integrator('DisplacementControl', 5, 1, -0.1)
ops.analysis('Static')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data = {}
data['step'] = []
data['load_factor'] = []
data['angle_degrees'] = []
print("Step,Load_Factor,Angle(degrees)")
M = 80
n_drawings = 10
for i in range(0, M+1):
    if i == 0:
        res = ops.analyze(0)
    else:
        res = ops.analyze(1)
    if res != 0:
        print(f"Analysis failed at step {i}")
        break
    lam = ops.getLoadFactor(1)
    theta = 2*np.pi-ops.eleResponse(hinges_element_numbers[-2], 'theta')[0]
    print(f"{i},{lam},{theta * 180.0 / np.pi}")
    data['step'].append(i)
    data['load_factor'].append(lam)
    data['angle_degrees'].append(theta * 180.0 / np.pi)

    if i % (M//n_drawings) == 0 or i == M:
        for node in range(len(vertices)):
            coor = ops.nodeCoord(node)
            disp = ops.nodeDisp(node)
            ax.scatter(coor[0]+disp[0], coor[1]+disp[1],
                       coor[2]+disp[2], color='b')
        for ele in ops.getEleTags():
            if ele in hinges_element_numbers:
                continue
            n1 = ops.eleNodes(ele)[0]
            n2 = ops.eleNodes(ele)[1]
            coor1 = ops.nodeCoord(n1)
            coor2 = ops.nodeCoord(n2)
            disp1 = ops.nodeDisp(n1)
            disp2 = ops.nodeDisp(n2)
            color = 'k'
            if i == M:
                color = 'r'
            ax.plot([coor1[0]+disp1[0], coor2[0]+disp2[0]],
                    [coor1[1]+disp1[1], coor2[1]+disp2[1]],
                    [coor1[2]+disp1[2], coor2[2]+disp2[2]], color=color)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_aspect('equal')

# Two legend lines, one red and one black
red_line = mlines.Line2D(
    [], [], color='r', label='Deformed shape at each step')
black_line = mlines.Line2D([], [], color='k', label='Elements')
plt.legend(handles=[red_line, black_line])
plt.show()

fig = plt.figure()
plt.plot(data['angle_degrees'], data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Angle (degrees)')
plt.ylabel('Load Factor')
plt.show()
a = 0
