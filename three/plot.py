#!/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data.csv")

for n in df["n"].unique():
    for type, ls in [("uni", "-"), ("cheba", "--")]:
        n_set = df[(df["n"] == n) & (df["type"] == type)]
        plt.plot(n_set["x"], n_set["err"], label=f"{type}, n={n}", linestyle=ls)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
