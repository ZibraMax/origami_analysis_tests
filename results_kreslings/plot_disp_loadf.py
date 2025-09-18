import numpy as np
import matplotlib.pyplot as plt
import json

filename = "results_kreslings/kresling_n_6_b_52_h_70_h0_60_t_0.5_kf_0.0.json"
with open(filename, 'r') as f:
    data = json.load(f)

solutions = data["solutions"]
lds = []
us = []
ts = []
for solution in solutions:
    info = solution['info']
    ld = round(info["ld"], 5)
    u = round(-info["vertical-disp"], 5)
    t = round(info["rotation_top_node"], 5)
    us.append(u)
    ts.append(t)
    lds.append(ld)
plt.plot(us, lds, 'r-', label='Vertical displacement')
plt.xlabel('Vertical displacement')
plt.ylabel('Load factor')
plt.grid()
fig = plt.figure()
plt.plot(ts, lds, 'b-', label='Rotation top node')
plt.ylabel('Load factor')
plt.xlabel('Rotation top node')
plt.grid()
plt.legend()
plt.show()
