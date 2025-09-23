import numpy as np
from unhinged import *
import opensees as ops
import matplotlib.pyplot as plt
import os


def model_setup():
    model = ShellAndHinge(O)
    model.add_material_shells(mat_tag=1, E=BASE_E, v=POSAO)
    model.add_material_shells(mat_tag=2, E=1e6*BASE_E, v=POSAO,
                              shell_list=list(range(2*n, 2*n+2)))
    model.add_material_bars(mat_tag=3, E=1e6*BASE_E, A=1.0)
    model.add_material_hinges(k=khinge)
    model.create_model()

    ops.fix(*[center_idx, 1, 1, 0, 0, 0, 0])

    model.setup_model(tol=TOL, maxiter=200)

    ops.integrator('DisplacementControl', center_idx, 3, -target_disp/M)
    # ops.integrator('MGDCM', 0.1, 15, 3, 1)
    ops.algorithm('Newton')
    ops.analysis('Static')
    return model


# Input parameters
H = 18.2
H0 = 0
n = 6
b = 13
thickness = 0.12
OPEN = True

BASE_E = 1219.4
POSAO = 0.3
khinge = 2.4e-3
mesh_refinement = 2
TOL = 1e-3
target_disp = H
M = 150
FOLDER = "results_otros"
EIGV_SCALE = 1
for EIGV_SCALE in np.linspace(0.0, 1, 11):
    try:
        os.mkdir(FOLDER)
    except Exception as e:
        pass

    Nmodes = 1
    number_sides = n
    basename = f"Kresling_n_{number_sides}_b_{b}_h_{H}_h0_{H0}_t_{thickness}_kf_{khinge}_poaso_{POSAO}_E_{BASE_E}_mrf_{mesh_refinement}_EIGVSCALE_{EIGV_SCALE}"
    if OPEN:
        basename = "Open" + basename
    kresling = Kresling(b=b, H0=H0, H1=H, n=n)
    data = kresling.generate(get_int_lines=False,
                             get_ext_lines=OPEN,
                             get_base_bars=False,
                             get_ext_hinges=False,
                             get_int_hinges=True,
                             get_panels=True,
                             get_base_panels=True,
                             get_base_hinges=True)
    delta_theta = kresling.delta_theta
    data['properties']['theta1'] = R(10)
    data['properties']['theta2'] = R(350)
    O = Geometry.from_json(data, t=thickness)
    O.mesh(n=mesh_refinement)

    res, center_idx = O.test_point_vertices([0.0, 0.0, H])
    O.add_bc_plane('z', 0.0, values=[1, 1, 1, 0, 0, 0])

    O.nbc.append([center_idx, 0, 0, -1, 0, 0, 0])

    model = model_setup()

    lam = ops.eigen('standard', 'symmBandLapack', Nmodes)
    eigenvectors = []
    for node in ops.getNodeTags():
        eigenvectors.append([])
        for mode in range(Nmodes):
            ev = ops.nodeEigenvector(node, mode+1)
            eigenvectors[-1].append(ev)

    model.solutions = []
    for mode in range(Nmodes):
        for node in ops.getNodeTags():
            nodedisp = ops.nodeDisp(node)
            for i, d in enumerate(nodedisp):
                ops.setNodeDisp(
                    node, i+1, eigenvectors[node][mode][i], '-commit')
        sol = model.get_disp_vector()
        sol["info"] = {"solver-type": "EIGEN", "ld": lam[mode]}
        model.solutions.append(sol)

    model.export_json(
        f"./{FOLDER}/eigv_{basename}.json")

    selected_solution = model.solutions[-1]['U']
    O.add_nodal_imperfections(selected_solution, scale=EIGV_SCALE)
    print(O.nbc)
    model = model_setup()

    # model.visualize()
    # plt.show()

    res = {"step": [], "load_factor": [], "disp": []}

    def callback(i):
        lam = ops.getLoadFactor(1)
        disp_z = ops.nodeDisp(center_idx, 3)
        rot_top = ops.nodeDisp(center_idx, 6)
        print(f"{i},{lam},{disp_z}")
        res['step'].append(i)
        res['load_factor'].append(lam)
        res['disp'].append(-disp_z)

        return {"solver-type": "DisplacementControl", "vertical-disp": disp_z, "rotation_top_node": rot_top}

    solutions = model.analyze(M, callback=callback)

    fig = plt.figure(figsize=[12, 6])
    ax2 = fig.add_subplot(1, 2, 1, projection='3d')
    ax = fig.add_subplot(1, 2, 2)
    model.visualize(ax=ax2)
    model.export_json(f"./{FOLDER}/{basename}.json")
    ax.plot(res['disp'], res['load_factor'], 'r-')
    ax.set_xlabel('Displacement')
    ax.set_ylabel('Load factor')
    ax.grid()
    plt.tight_layout()
    # plt.savefig(f"./{FOLDER}/figures/{basename}.pdf", dpi=300)
    plt.show()
