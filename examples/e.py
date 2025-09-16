import unhinged as uh
import numpy as np
import matplotlib.pyplot as plt


# Create hexagon coordinates
def hexagon(R):
    nodes = []
    coords = []
    for i in range(6):
        angle = np.pi/3 * i
        x = R * np.cos(angle)
        y = R * np.sin(angle)
        coords.append([x, y, 0.0])
        nodes.append(i)

    return coords, nodes


O = uh.Geometry(6)
coords, nodes = hexagon(1.0)
O.add_nodes(coords)
hexagon_panel = uh.PolygonalPanel(nodes, thickness=0.1)
O.add_panel(hexagon_panel)
O.mesh(n=2)
model = uh.ShellAndHinge(O)
model.add_material_shells(mat_tag=1, E=200, v=0.3)
model.add_material_hinges(k=0.5)
model.create_model()
model.setup_model(tol=1e-2)
model.visualize(undeformed=False, node_labels=True)
plt.show()
