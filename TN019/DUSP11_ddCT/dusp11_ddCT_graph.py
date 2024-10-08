import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


infile = pd.read_excel("ddCT_DUSP11_table_1.xlsx")

sns.set(style="ticks")
plot = sns.barplot(
    data=infile,
    x="Sample Name",
    y="value",
    errorbar="sd",
    hue = "Sample Name",
    palette="Paired",
    legend=False
)
plt.xlabel(" ")
plt.ylabel("foldchange")
plt.ylim(0,1)
plt.legend(title="")
plt.title("DUSP11 mRNA foldchange during infection")
plt.savefig("TN019.21_DUSP11_ddCT.png", dpi=300, format="png")
plt.savefig("TN019.21_DUSP11_ddCT.svg", dpi=300, format="svg")