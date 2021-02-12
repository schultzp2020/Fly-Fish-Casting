import pandas as pd
import numpy as np


t = np.linspace(0., 0.1, 11)

xt = np.zeros(len(t))
# xt_dict = {round(t[i], 2): xt[i] for i in range(len(t))}
yt = 0.5*np.ones(len(t))

xt_dot = np.zeros(len(t))
yt_dot = np.zeros(len(t))

xt_ddot = np.zeros(len(t))
yt_ddot = np.zeros(len(t))

data = {"Time Intervals": t, "xt": xt, "yt": yt, "xt dot": xt_dot, "yt dot": yt_dot, "xt ddot": xt_ddot, "yt ddot": yt_ddot}
columns = ["Time Intervals", "xt", "yt", "xt dot", "yt dot", "xt ddot", "yt ddot"]

df = pd.DataFrame(data, columns=columns)

df.to_csv("Flyrod.csv", index=False)

print(df)