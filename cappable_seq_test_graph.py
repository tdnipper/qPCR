import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Import data, dropna and remove NTC
data = pd.read_excel("cappable_seq_test_output.xlsx").dropna()
data = data.loc[~data['Sample Name'].str.contains('NTC')]

# Create catplot separating by Sample Name and showing each target
# g = sns.catplot(
#     data=data,
#     kind="bar",
#     x="Target Name",
#     y="percent_recovery",
#     hue="Sample Name",
#     col="Sample Name",
#     col_wrap=4,
#     height=4,
#     aspect=1,
#     palette="muted",
# )

# T-test between mock and ix for each Target Name
grouped_data = data.groupby(['Sample Name', 'Target Name'])['percent_recovery']

g = sns.barplot(
    data=data, x="Target Name", y="percent_recovery", hue="Sample Name", palette="muted"
)

# Rotate x-axis labels for all plots
# Adjust the y-axis to 25
# Adjust the margins
plt.gca().set_xticklabels(data["Target Name"].unique(), rotation=45, ha="right")
g.set(ylim=(0, 25))
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9)


# Save the plot
sns.set_style("whitegrid")
plt.savefig("cappable_seq_test.png", format="png", dpi=300)
