import numpy as np
import matplotlib.pyplot as plt


class Kresling():
    def __init__(self, b, H0, H1, n):
        self.b = b
        self.H0 = H0
        self.H1 = H1
        self.n = n
        pi_n = np.pi / n
        cot_pi_n = 1 / np.tan(pi_n)
        if abs(H1**1-H0**2) <= b**2*cot_pi_n**2:
            Exception(
                "The pattern is not valid for the given parameters. Contact between the panels")

        h0 = H0 / b
        h = H1 / b

        sin_pi_n = np.sin(pi_n)
        cos_pi_n = np.cos(pi_n)
        csc_pi_n = 1 / sin_pi_n
        numerator = 2*sin_pi_n * \
            (sin_pi_n*(cot_pi_n**2*csc_pi_n**2-(h**2-h0**2)**2)**0.5-cos_pi_n)
        x1 = numerator/(1+h**1-h0**2+(1-h**2+h0**2)*np.cos(2*np.pi/n))
        x2 = numerator/(1-h**1+h0**2+(1+h**2-h0**2)*np.cos(2*np.pi/n))
        alpha = np.arccos((x2*(x2-cot_pi_n)) /
                          (((x2**2+1)*(h0**2*(x2**2+1)+x2**2*csc_pi_n**2))**0.5))

        phi1 = 2*np.arctan(x1)
        phi0 = 2*np.arctan(x2)
        self.delta_theta = phi1 - phi0

        self.props = {
            "phi0": phi0,
            "phi1": phi1,
            "alpha": alpha
        }

    def rotate_z_axis(self, point, angle):
        """Rotate a point around the Z-axis by a given angle in radians."""
        x, y, z = point
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        x_rotated = x * cos_angle - y * sin_angle
        y_rotated = x * sin_angle + y * cos_angle
        return [float(x_rotated), float(y_rotated), float(z)]

    def generate(self, get_panels=True, get_int_hinges=True, get_ext_hinges=True, get_base_bars=True, get_int_lines=True, get_ext_lines=True, get_base_panels=True, get_base_hinges=True):
        design = self.props
        phi0 = design["phi0"]
        phi1 = design["phi1"]
        alpha = design["alpha"]
        b = self.b
        H0 = self.H0
        H1 = self.H1
        n = self.n

        kresling = {}
        kresling["nodes"] = []
        kresling["triangles"] = []
        kresling["base_lines"] = []
        kresling["interior_lines"] = []
        kresling["exterior_lines"] = []
        kresling["interior_hinges"] = []
        kresling["exterior_hinges"] = []
        kresling["panel_hinges"] = []
        kresling["base_panels"] = []

        nodes = kresling["nodes"]
        triangles = kresling["triangles"]
        ext_lines = kresling["exterior_lines"]
        int_lines = kresling["interior_lines"]
        int_hinges = kresling["interior_hinges"]
        ext_hinges = kresling["exterior_hinges"]
        base_lines = kresling["base_lines"]
        panel_hinges = kresling["panel_hinges"]
        base_panels = kresling["base_panels"]

        # Create base and top polygons in 3D coordinates
        for i in range(n):
            angle = 2 * np.pi * i / n
            x_base = b * np.cos(angle)
            y_base = b * np.sin(angle)
            nodes.append([float(x_base), float(y_base), 0.0])

        for i in range(n):
            angle = 2 * np.pi * i / n
            x_base = b * np.cos(angle)
            y_base = b * np.sin(angle)
            x_top = x_base
            y_top = y_base
            coodr = (x_top, y_top, H1)
            coodr = self.rotate_z_axis(coodr, phi1)
            nodes.append(coodr)

        m = len(nodes)
        for i in range(n-1):
            xi = i
            xip1 = (i + 1)
            xipn = (n + i)
            xipn1 = (n + (i + 1))
            triangles.append([xi, xip1, xipn1])
            triangles.append([xi, xipn1, xipn])

        for i in range(n):
            base_lines.append([i, (i + 1) % n])           # Base edges
        for i in range(n):
            base_lines.append([n + i, n + (i + 1) % n])   # Top edges

        # add panel lines
        for i in range(n):
            ext_lines.append([i, n + i])
            int_lines.append([i, n + (i + 1) % n])

        for i in range(n-1):
            I = (i-1) % n
            J = (i) % n
            K = (i+n) % (2*n)
            L = (i+n+1) % (2*n)
            ext_hinges.append([I, J, K, L])
        ext_hinges.append([n-2, n-1, 2*n-1, n])

        for i in range(n-1):
            I = (i+n) % (2*n)
            J = (i+1+n) % (2*n)
            K = (i) % (n)
            L = (i+1) % (n)
            int_hinges.append([I, J, K, L])
        int_hinges.append([2*n-1, n, n-1, 0])

        triangles.append([n-1, 0, n])
        triangles.append([n-1, n, m-1])

        # m = len(nodes)+1
        # for i in range(n):
        #     panel_hinges.append([m-1, i, (i+1) % n, (i+n+1) % (2*n)])

        # m = len(nodes)+2
        # for i in range(n):
        #     panel_hinges.append([m-1, (i+n) % (2*n), (i+n+1) % (2*n), i])

        # kresling["dictionary"] = triangles + ext_lines + \
        #     int_lines + int_hinges + ext_hinges + base_lines
        # kresling["types"] = ["Shell"]*len(triangles) + \
        #     ["Opening"]*len(ext_lines) + \
        #     ["Bar"]*len(int_lines) + \
        #     ["OriHinge"]*len(int_hinges) + \
        #     ["OriHinge"]*len(ext_hinges) + ["Bar"]*len(base_lines)

        kresling["dictionary"] = []
        kresling["types"] = []
        if get_panels:
            kresling["dictionary"] += triangles
            kresling["types"] += ["Shell"]*len(triangles)
        if get_ext_lines:
            kresling["dictionary"] += ext_lines
            kresling["types"] += ["Opening"]*len(ext_lines)
        if get_int_lines:
            kresling["dictionary"] += int_lines
            kresling["types"] += ["Bar"]*len(int_lines)
        if get_int_hinges:
            kresling["dictionary"] += int_hinges
            kresling["types"] += ["OriHinge"]*len(int_hinges)
        if get_ext_hinges:
            kresling["dictionary"] += ext_hinges
            kresling["types"] += ["OriHinge"]*len(ext_hinges)
        if get_base_bars:
            kresling["dictionary"] += base_lines
            kresling["types"] += ["Bar"]*len(base_lines)

        if get_base_panels:

            base1 = [i for i in range(n)]
            base2 = [i for i in range(n, 2*n)]
            base_panels.append(base1)
            base_panels.append(base2)

            kresling["dictionary"] += base_panels
            kresling["types"] += ["Poly"]*len(base_panels)

        kresling["properties"] = {"problem": "ShellAndHinge", }

        return kresling

    def test_point_vertices(self, vertex, vertices, tol=1e-10):
        if len(vertices) == 0:
            return False, None
        vertex = np.array(vertex)
        vertices = np.array(vertices)

        distances = np.linalg.norm(vertices-vertex, axis=1)
        if np.min(distances) < tol:
            return True, int(np.argmin(distances))
        return False, None
