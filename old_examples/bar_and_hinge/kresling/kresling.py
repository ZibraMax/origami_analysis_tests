import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib
matplotlib.use('TkAgg')

A_MODIFIER_STIFF = 1e8
ELE_TYPES = {}
LOADS = {}


def get_disp_vector():
    U_NODES = []
    nodes = ops.getNodeTags()
    lam = ops.getLoadFactor(1)
    for j in nodes:
        U_NODES.extend(ops.nodeDisp(j))
    U_NODES = np.array(U_NODES)
    U_NODES = U_NODES.reshape((3*len(nodes), 1))
    # Get integrator name
    return {
        "info": {"solver-type": "int_name", "ld": lam},
        "U": U_NODES.tolist(),
    }


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
        # Plot loads as unit lengt arrows in plot coords
        for node, load in LOADS.items():
            coor = np.array(ops.nodeCoord(node))
            load = np.array(load)
            if np.linalg.norm(load) > 0:
                load = load / np.linalg.norm(load) * 0.1 * np.linalg.norm(np.array(ops.nodeCoord(
                    ops.getNodeTags()[-1])) - np.array(ops.nodeCoord(ops.getNodeTags()[0])))
                coor = coor - load
                ax.quiver(coor[0], coor[1], coor[2], load[0], load[1],
                          load[2], color='red', arrow_length_ratio=0.5, linewidth=2)

        ax.axis('off')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')


def create_model_from_json(file_path):
    data = json.load(open(file_path, 'r'))
    nodes = np.array(data['nodes'])
    dictionary = data['dictionary']
    types = data['types']
    interior_bars = data['interior_bars']
    exterior_bars = data['exterior_bars']
    base_bars = data['base_bars']
    ebc = data['ebc']
    nbc = data['nbc']
    theta_1 = data['theta_1']
    theta_2 = data['theta_2']
    nvn = data['nvn']
    properties = data['properties']

    E = properties['E']
    A = properties['A']
    kf = properties['kf']

    if not isinstance(A, list):
        A = [A]*len(dictionary)
    if not isinstance(E, list):
        E = [E]*len(dictionary)
    if not isinstance(kf, list):
        kf = [kf]*len(dictionary)
    assert len(A) == len(
        dictionary), f"Length of A ({len(A)}) does not match number of elements ({len(dictionary)})"
    assert len(E) == len(
        dictionary), f"Length of E ({len(E)}) does not match number of elements ({len(dictionary)})"
    assert len(kf) == len(
        dictionary), f"Length of kf({len(kf)}) does not match number of elements ({len(dictionary)})"

    data["properties"]["E"] = E
    data["properties"]["kf"] = kf
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', nvn)

    materials = np.unique(E).tolist()
    for i, e in enumerate(materials):
        ops.uniaxialMaterial("Elastic", i, float(e))

    element_to_material = {}
    for i, e in enumerate(E):
        idx = materials.index(e)
        element_to_material[i] = idx

    for i, node in enumerate(nodes):
        ops.node(i, *node)
    for i, (elem_nodes, etype) in enumerate(zip(dictionary, types)):
        nel = len(ops.getEleTags())
        ELE_TYPES[nel] = etype
        if etype == "OriHinge":
            ops.element("OriHinge", nel, *elem_nodes, kf[i], theta_1, theta_2)
        elif etype == "CorotTruss" or etype == "Truss":
            if elem_nodes in base_bars:
                A[i] = A[i]*A_MODIFIER_STIFF
            ops.element(etype, nel, *elem_nodes, A[i], element_to_material[i])
    data["properties"]["A"] = A

    for bc in ebc:
        ops.fix(*bc)

    return data, materials


if __name__ == '__main__':
    file_path = "./bar_and_hinge/kresling/kresling.json"
    data, materials = create_model_from_json(file_path)
    data["solutions"] = []
    base_bars = data['base_bars']
    nodes_base = []
    for bar in base_bars:
        nodes_base.extend(bar)
    nodes_base = list(set(nodes_base))
    nodes_bottom_base = []
    for node in nodes_base:
        if data['nodes'][node][2] == 0:
            nodes_bottom_base.append(node)
    nodes_top_base = []
    for node in nodes_base:
        if data['nodes'][node][2] == max(np.array(data['nodes'])[:, 2]):
            nodes_top_base.append(node)

    # Fix bottom base nodes
    for node in nodes_bottom_base:
        ops.fix(node, 1, 1, 1)

    # Equal Z dof for top base nodes
    for i in range(1, len(nodes_top_base)):
        ops.equalDOF(nodes_top_base[0], nodes_top_base[i], 3)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    P = 1
    for node in nodes_top_base:
        ops.load(node, 0, 0, -P/len(nodes_top_base))
        LOADS[node] = [0, 0, -P/len(nodes_top_base)]

    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1.0e-6, 100)

    M = 100
    ops.integrator('DisplacementControl', nodes_top_base[0], 3, -20/M)
    # ops.integrator('LoadControl', 0.7)
    ops.integrator("MGDCM", 0.1, 6, 2, 1)
    ops.algorithm('Newton')
    ops.analysis('Static')

    fig = plt.figure(figsize=[6, 6])
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    visualize(ax=ax, undeformed=True, node_labels=True)
    plt.show()

    fig = plt.figure(figsize=[12, 6])
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    N = 1
    nodes_with_loads = list(LOADS.keys())
    res = {"step": [], "load_factor": [], "disp": []}
    try:
        for i in range(0, M+1):
            if i != 0:
                if ops.analyze(1) != 0:
                    print(f"Analysis failed at step {i}")
                    break
            lam = ops.getLoadFactor(1)
            disp_z = ops.nodeDisp(nodes_with_loads[0], 3)
            print(f"{i},{lam},{disp_z}")
            res['step'].append(i)
            res['load_factor'].append(lam)
            res['disp'].append(-disp_z)

            if i % (M//N) == 0:
                visualize(ax=ax, plot_hinges=False)
            data["solutions"].append(get_disp_vector())
    except KeyboardInterrupt as e:
        print(f"Analysis stopped at step {i} with error: {e}")

    # Change the types as follows: CorotTruss -> L1V, OriHinge -> OH

    data["types"] = ["L1V" if t == "CorotTruss" else "OH" for t in data["types"]]

    json.dump(data, open('./output/kresling_bah.json', 'w'))
    visualize(ax=ax, plot_hinges=False)

    ax2.plot(res['disp'], res['load_factor'], 'r-')
    ax2.set_xlabel('Displacement')
    ax2.set_ylabel('Load factor')
    ax2.grid()

    plt.tight_layout()
    plt.show()
