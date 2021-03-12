import numpy as np
import pandas as pd

def get_y(phi, i, yt, length):
    #get y at position i
    #rod tip y coordinate at yt
    #L is vector of segment Lengths
    s = yt
    for j in range(i):
        s += length[j] * np.sin(phi[j])
    s+= length[i] * np.sin(phi[i]) * 0.5 
    return s

npts = 60

length = 0.1 * np.ones(npts)
radius = 0.0005205 * np.ones(npts)
mass = 1158*np.pi* np.power(0.0005205,2) * 0.1 * np.ones(npts)
phi_init = (-np.pi/2) * np.ones(npts)
phi_dot_init = np.zeros(npts)

for i in range(npts):
  y_pos = get_y(phi_init, i, 3, length)
  if y_pos < 0:
    phi_init[i] = 0.0

data = {"Length": length, "Radius": radius, "Mass": mass, "Phi Init": phi_init, "Phi Dot Init": phi_dot_init}
columns = ["Length", "Radius", "Mass", "Phi Init", "Phi Dot Init"]

df = pd.DataFrame(data, columns=columns)

df.to_csv("Flyline.csv", index=False)

print(df)