#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_excel("TN042.2.xlsx")

dot = ["#87CEEB"]
bar = ["#4682B4"]
g = sns.swarmplot(
    x="Virus",
    y="Titer",
    data=data,
    hue="Virus",
    palette=dot,
    edgecolor="black",
    linewidth=0.5,
)
g = sns.barplot(
    x="Virus", y="Titer", data=data, hue="Virus", errorbar=None, palette=bar
)
plt.yscale("log")
plt.ylim(1e5, 1e7)
plt.title("WSN Titer during dox induction")

plt.savefig("TN042.2.png", dpi=300)
