import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib
matplotlib.use('TkAgg')
EXTRA_NODES = True
A_MODIFIER_STIFF = 1e10
ELE_TYPES = {}
LOADS = {}
t = 10


def test_point_vertices(vertex, vertices, tol=1e-10):
    if len(vertices) == 0:
        return False, None
    vertex = np.array(vertex)
    vertices = np.array(vertices)

    distances = np.linalg.norm(vertices-vertex, axis=1)
    if np.min(distances) < tol:
        return True, np.argmin(distances)
    return False, None


def get_disp_vector():
    U_NODES = []
    nodes = ops.getNodeTags()
    lam = ops.getLoadFactor(1)
    for j in nodes:
        U_NODES.extend(ops.nodeDisp(j))
    U_NODES = np.array(U_NODES)
    U_NODES = U_NODES.reshape((6*len(nodes), 1))
    # Get integrator name
    return {
        "info": {"solver-type": "int_name", "ld": lam},
        "U": U_NODES.tolist(),
    }


def create_extra_nodes_for_hinges(nodes, elemements, types):
    if not isinstance(nodes, np.ndarray):
        nodes = np.array(nodes)
    nodes = nodes.tolist()
    nodes_hinges_middle = []
    nodes_to_element = {}
    tie_nodes = []
    for i, etype in enumerate(types):
        if etype == "OriHinge":
            n1, n2, n3, n4 = elemements[i]
            nodes_hinges_middle.append(n2)
            nodes_hinges_middle.append(n3)
            # Center of n2 and n3
            new_node = (np.array(nodes[n2]) + np.array(nodes[n3])) / 2
            res, idx_new_node = test_point_vertices(new_node, nodes, tol=1e-2)
            if res:
                nodes_hinges_middle.append(int(idx_new_node))
            else:
                print(new_node, idx_new_node)
                raise ValueError("Hinge middle node not found in nodes list")

        else:
            element = elemements[i]
            for n in element:
                if n not in nodes_to_element:
                    nodes_to_element[n] = []
                nodes_to_element[n].append(i)
    nodes_hinges_middle = list(set(nodes_hinges_middle))
    for node in nodes_hinges_middle:
        connected_elements = nodes_to_element[node]
        number_conected = len(connected_elements)
        for i in range(number_conected - 1):
            new_node_id = len(nodes)
            nodes.append(nodes[node])
            tie_nodes.append((node, new_node_id))
            elem_idx = connected_elements[i+1]
            pos_new_node = elemements[elem_idx].index(node)
            elemements[elem_idx][pos_new_node] = new_node_id

    return np.array(nodes), elemements, tie_nodes


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

    if EXTRA_NODES:
        nodes, dictionary, tie_nodes = create_extra_nodes_for_hinges(
            nodes, dictionary, types)
        data['nodes'] = nodes.tolist()
        data['dictionary'] = dictionary
        data['tie_nodes'] = tie_nodes
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
    v = properties['v']
    kf = properties['kf']

    if not isinstance(v, list):
        v = [v]*len(dictionary)
    if not isinstance(E, list):
        E = [E]*len(dictionary)
    if not isinstance(kf, list):
        kf = [kf]*len(dictionary)
    assert len(v) == len(
        dictionary), f"Length of v ({len(v)}) does not match number of elements ({len(dictionary)})"
    assert len(E) == len(
        dictionary), f"Length of E ({len(E)}) does not match number of elements ({len(dictionary)})"
    assert len(kf) == len(
        dictionary), f"Length of kf({len(kf)}) does not match number of elements ({len(dictionary)})"

    data["properties"]["E"] = E
    data["properties"]["v"] = v
    data["properties"]["kf"] = kf
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', nvn)

    materials = list(set(list(zip(E, v))))
    for i, e in enumerate(materials):
        ops.section('ElasticMembranePlateSection', i, *e, t)

    ops.uniaxialMaterial("Elastic", 999, 1)
    element_to_material = {}
    for i, e in enumerate(zip(E, v)):
        idx = materials.index(e)
        element_to_material[i] = idx

    for i, node in enumerate(nodes):
        ops.node(i, *node)
    for i, (elem_nodes, etype) in enumerate(zip(dictionary, types)):
        nel = len(ops.getEleTags())
        ELE_TYPES[nel] = etype
        if etype == "OriHinge":
            ops.element("OriHinge", nel, *elem_nodes, kf[i], theta_1, theta_2)
        elif etype == "ShellT" or etype == "ASDShellT3":
            ops.element("ASDShellT3", nel, *elem_nodes,
                        element_to_material[i], '-corotational')
    for i, e in enumerate(base_bars):
        # Add bars as infinite stiff
        nel = len(ops.getEleTags())
        ELE_TYPES[nel] = "CorotTruss"
        ops.element("CorotTruss", nel, *e, A_MODIFIER_STIFF, 999)

    for bc in ebc:
        ops.fix(*bc)
    if EXTRA_NODES:
        for bc in tie_nodes:
            if bc[0] != bc[1]:
                print(bc)
                ops.equalDOF(bc[0], bc[1], 1, 2, 3)

    return data, materials


if __name__ == '__main__':
    file_path = "./shell_and_hinge/kresling/kresling.json"
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
        ops.fix(node, 1, 1, 1, 0, 0, 0)

    # Equal Z dof for top base nodes
    # for i in range(1, len(nodes_top_base)):
    #     ops.equalDOF(nodes_top_base[0], nodes_top_base[i], 3)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    P = 1
    for node in nodes_top_base:
        ops.load(node, 0, 0, -P/len(nodes_top_base), 0, 0, 0)
        LOADS[node] = [0, 0, -P/len(nodes_top_base)]

    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1.0e-6, 100)

    M = 250
    # ops.integrator('DisplacementControl', nodes_top_base[0], 3, -20/M)
    # ops.integrator('LoadControl', 0.7)
    ops.integrator("MGDCM", 10, 6, 3, 1)
    ops.algorithm('Newton')
    ops.analysis('Static')

    fig = plt.figure(figsize=[6, 6])
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    visualize(ax=ax, undeformed=True, node_labels=True)
    plt.show()

    fig = plt.figure(figsize=[12, 6])
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    N = 2
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

    data["types"] = ["T1V" if t == "ASDShellT3" else "OH" for t in data["types"]]

    json.dump(data, open('./output/rigid_kresling_sah.json', 'w'))
    visualize(ax=ax, plot_hinges=False)

    ax2.plot(res['disp'], res['load_factor'], 'r-')
    ax2.set_xlabel('Displacement')
    ax2.set_ylabel('Load factor')
    ax2.grid()

    plt.tight_layout()
    plt.show()
