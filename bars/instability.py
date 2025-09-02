"""Example of instability (truss) non linear geometric with snap-through. usign opensess"""
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib
datas = []

for integrator in ["DisplacementControl", "MGDCM"]:
    print(f"Using integrator: {integrator}")
    # Geometry and material
    a = 400
    b = 20
    L = (a**2 + b**2)**0.5
    h = 1
    t = 653
    A = h*t
    P = 1
    young = 20500  # Young's modulus in MPa
    EA = young*A

    coords = [
        [0, 0],
        [a, b],
        [2*a, 0]
    ]
    elements = [[0, 1],
                [1, 2]]

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.uniaxialMaterial('Elastic', 1, young)

    for node_id, coord in enumerate(coords):
        ops.node(node_id, *coord)

    for elem_id, (n1, n2) in enumerate(elements):
        ops.element('CorotTruss', elem_id, n1, n2, A, 1)

    ops.fix(0, 1, 1)
    ops.fix(2, 1, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(1, 0, -P)

    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1.0e-3, 100)
    # ops.integrator('LoadControl', 0.1)
    if integrator == "MGDCM":
        ops.integrator('MGDCM', 80)
    else:
        ops.integrator('DisplacementControl', 1, 2, -0.5)
    ops.algorithm('Newton')
    ops.analysis('Static')

    data = {}
    data['step'] = []
    data['load_factor'] = []
    data['disp'] = []

    M = 90

    try:
        for i in range(0, M+1):
            if i != 0:
                if ops.analyze(1) != 0:
                    print(f"Analysis failed at step {i}")
                    break
            else:
                print(f"Analysis started")
            lam = ops.getLoadFactor(1)
            disp_z = ops.nodeDisp(1, 2)
            print(f"{i},{lam},{disp_z}")
            data['step'].append(i)
            data['load_factor'].append(lam)
            data['disp'].append(-disp_z)
    except KeyboardInterrupt as e:
        print(f"Analysis stopped at step {i} with error: {e}")
    datas.append((integrator, data))

for integrator, data in datas:
    plt.plot(data['disp'], data['load_factor'], 'o--', label=integrator)

plt.xlabel('Mid node vertical displacement')
plt.ylabel('Load factor')
plt.legend()
plt.grid()
plt.title('Snap-through instability of a truss')
plt.show()
