import matplotlib.lines as mlines
import opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import json


def R(degrees): return degrees * np.pi / 180.0
def D(radians): return radians * 180.0 / np.pi


THETA_1 = 10
THETA_2 = 350

data = json.load(open('./shell_and_hinge/simple_fold/fold.json'))

nodes = data['nodes']
dictionary = data['dictionary']
types = data['types']
ebc = data['ebc']
nbc = data['nbc']
nvn = data['nvn']
props = data['properties']
E = props['E']
v = props['v']
t = props['t']
kf = props['kf']


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
ops.section('ElasticMembranePlateSection', 1, E, v, t)
