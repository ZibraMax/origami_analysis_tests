from unhinged import *

H = 70
n = 6
b = 68

data = Kresling(b=b, H0=0.0, H1=H, n=n).generate()
print(data)
O = Geometry.from_json(data, t=0.2)
