import matplotlib.lines as mlines
import opensees as ops
import math
import matplotlib.pyplot as plt


def R(degrees): return degrees * math.pi / 180.0
def D(radians): return radians * 180.0 / math.pi


# --- Test OriHinge ---
ops.wipe()
ops.model('BasicBuilder', '-ndm', 3, '-ndf', 6)

# --- Parameters ---
L = 1.0
theta_0 = 210
phi = 2.0 * math.pi - theta_0 * math.pi / 180.0
x_coord = math.sin(math.pi / 3.0) * L
z_coord = math.sin(phi) * x_coord
x_coord2 = math.cos(phi) * x_coord

# --- Nodes ---
ops.node(0, 0.0, -0.5*L, 0.0)
ops.node(1, 0.0, 0.5*L, 0.0)
ops.node(2, x_coord, 0.0, 0.0)
ops.node(3, x_coord2, 0.0, z_coord)

# --- Supports ---
ops.fix(0, 1, 1, 1, 1, 1, 1)
ops.fix(1, 1, 1, 1, 1, 1, 1)
ops.fix(2, 1, 1, 1, 1, 1, 1)
ops.fix(3, 0, 0, 0, 1, 1, 1)

# --- Material ---
ops.uniaxialMaterial('Elastic', 1, 1000.0)

# --- Elements ---
ops.element('corotTruss', 0, 0, 1, 1.0, 1)
ops.element('corotTruss', 1, 1, 2, 1.0, 1)
ops.element('corotTruss', 2, 2, 0, 1.0, 1)
ops.element('corotTruss', 3, 0, 3, 1.0, 1)
ops.element('corotTruss', 4, 1, 3, 1.0, 1)
ops.element('OriHinge', 5, 3, 0, 1, 2, 0.3, R(10), R(350))

# --- Nodal Load ---
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# --- Static Analysis Setup ---
ops.system('Umfpack')
ops.constraints('Plain')
ops.numberer('RCM')
ops.test('NormDispIncr', 1.0e-3, 200)
# ops.algorithm('Newton')
ops.integrator('ArcLength', 0.1, 0.1)
ops.analysis('Static')

# --- Perform Analysis ---
data = {}
data['step'] = []
data['load_factor'] = []
data['angle_degrees'] = []
data['disp_node4_z'] = []
print("Step,Load_Factor,Angle(degrees),Disp_Node4_Z")
fig = plt.figure(figsize=[12, 5])
ax = fig.add_subplot(1, 2, 1, projection='3d')
M = 20
hinges_element_numbers = [5]
for i in range(1, M+1):
    if ops.analyze(1) != 0:
        print(f"Analysis failed at step {i}")
        break
    lam = ops.getLoadFactor(1)
    theta = ops.eleResponse(5, 'theta')
    disp_z = ops.nodeDisp(3, 3)
    print(f"{i},{lam},{theta[0] * 180.0 / math.pi},{disp_z}")
    data['step'].append(i)
    data['load_factor'].append(lam)
    data['angle_degrees'].append(theta[0] * 180.0 / math.pi)
    data['disp_node4_z'].append(disp_z)
    for node in range(4):
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
    [], [], color='r', label='Final shape')
black_line = mlines.Line2D([], [], color='k', label='Deformation process')
plt.legend(handles=[red_line, black_line])
ax = fig.add_subplot(1, 2, 2)
plt.plot(data['angle_degrees'], data['load_factor'], 'o--')
plt.grid()
plt.xlabel('Angle (degrees)')
plt.ylabel('Load Factor')
plt.tight_layout()
plt.show()
