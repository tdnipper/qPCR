import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_excel("TN022.1_output.xlsx")
plot = sns.FacetGrid(data, col="Sample Name")
plot.map(
    sns.histplot,
    data=data,
    x="Target Name",
    y="percent_recovery",
)
plot.set_xticklabels(rotation=45)
plt.subplots_adjust(bottom=0.3)
plt.show()
