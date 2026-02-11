import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_excel("TN022.2_output.xlsx")

# Assign group name based on timepoint and MOI
data["group"] = (
    data["Sample Name"]
    .astype(str)
    .apply(
        lambda x: "MOI5 T12"
        if x.startswith("MOI5 T12")
        else "MOI5 T16"
        if x.startswith("MOI5 T16")
        else "MOI5 T24"
        if x.startswith("MOI5 T24")
        else "MOI0.02 T24"
    )
)

# Assign infection status based on mock or ix in sample name
data["infect"] = (
    data["Sample Name"]
    .astype(str)
    .apply(lambda x: "mock" if x.endswith("mock") else ("ix"))
)
colors = ["#00a9e0", "#00549f"]

# Create a faceted plot separating by group name and coloring based on infection status
g = sns.catplot(
    data=data,
    kind="bar",
    x="Target Name",
    y="percent_recovery",
    hue="infect",
    col="group",
    palette=colors,
    col_wrap=4,
    height=4,
    aspect=1,
    legend=False,
)

# Show the plot
sns.despine()
sns.set_style("whitegrid")

# plt.xticks(rotation=45)
g.set(ylim=(0, 15))
plt.savefig("TN022.2_sep_output.png", format="png", dpi=300)
