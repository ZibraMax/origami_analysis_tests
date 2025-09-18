import numpy as np
import matplotlib.pyplot as plt
import json
import os
import pandas as pd
FOLDER = "./results_kreslings/"
files = os.listdir(FOLDER)
filenames = []
ns = []
bs = []
hs = []
h0s = []
ts = []
kfs = []
for file in files:
    if file.endswith(".json") and not file.startswith("eigv"):
        # exmaple_filename = kresling_n_8_b_52_h_70_h0_30_t_0.5_kf_0.0
        n = int(file.split("_")[2])
        b = float(file.split("_")[4])
        h = float(file.split("_")[6])
        h0 = float(file.split("_")[8])
        t = float(file.split("_")[10])
        kf = float(file.split("_")[12].replace(".json", ""))
        filenames.append(f"{FOLDER}/{file}")
        ns.append(n)
        bs.append(b)
        hs.append(h)
        h0s.append(h0)
        ts.append(t)
        kfs.append(kf)
df = pd.DataFrame({"filename": filenames, "n": ns, "b": bs,
                  "h": hs, "h0": h0s, "t": ts, "kf": kfs})

df1 = df[df.n == 8]

df1.sort_values(by="kf", inplace=True)

for filename, kf in zip(df1["filename"], df1["kf"]):
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
    plt.plot(us, lds, '-', label=f'k={kf}')
plt.xlabel('Vertical displacement')
plt.ylabel('Load factor')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
