import numpy as np
import pandas as pd

npts = 10

length = 0.1 * np.ones(npts)
radius = 0.001 * np.ones(npts)
mass = 0.01 * np.ones(npts)
phi_init = np.zeros(npts)
phi_dot_init = np.zeros(npts)

data = {"Length": length, "Radius": radius, "Mass": mass, "Phi Init": phi_init, "Phi Dot Init": phi_dot_init}
columns = ["Length", "Radius", "Mass", "Phi Init", "Phi Dot Init"]

df = pd.DataFrame(data, columns=columns)

df.to_csv("Flyline.csv", index=False)

print(df)