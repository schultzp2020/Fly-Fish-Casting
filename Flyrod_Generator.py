import pandas as pd
import numpy as np

t = np.linspace(0., 1.0, 101)

data = {"Time Intervals": t}
columns = ["Time Intervals"]

df = pd.DataFrame(data, columns=columns)

df.to_csv("Flyrod.csv", index=False)

print(df)