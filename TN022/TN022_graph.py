import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_excel("TN022.1_output.xlsx")

# # Create barplot using facet grid separated by Sample Name
# plot = sns.FacetGrid(data,
#                      col="Sample Name",
#                      height=4,
#                      aspect=2)
# plot.map(sns.barplot,
#          data=data,
#          x="Target Name",
#          y="percent_recovery")
# plot.set_xticklabels(rotation=45)
# plt.subplots_adjust(bottom=0.3)
# plt.show()

# Create catplot separating by Sample Name and showing each target
g = sns.catplot(
    data=data,
    kind="bar",
    x="Target Name",
    y="percent_recovery",
    hue="Sample Name",
    col="Sample Name",
    col_wrap=2,
    height=4,
    aspect=1,
    palette="muted",
)

# Rotate x-axis labels for all plots
g.set_xticklabels(rotation=45, horizontalalignment="right")

# Show the plot
sns.despine()
sns.set_style("whitegrid")
# plt.show()
plt.savefig("TN022.1_sep_output.png", format="png", dpi=300)
