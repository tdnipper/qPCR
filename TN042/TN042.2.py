#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg

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

# Calculate the p-value
ttest = pg.ttest(
    data["Titer"][data["Virus"] == "CRISPRi - Dox"],
    data["Titer"][data["Virus"] == "CRISPRi + Dox"],
    paired=False,
    alternative="two-sided",
)
ttest.to_excel("TN042.2_ttest.xlsx")
def annotate_stats(pvalue, y, h):
    plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c="k")
    plt.text(0.5, y + 2 * h, f"p = {pvalue:.2}", ha="center", va="center")
annotate_stats(ttest["p-val"].values[0], 5e6, 5e5)

plt.savefig("TN042.2.png", dpi=300)
