import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# Import data, dropna and remove NTC
data = pd.read_excel("cappable_seq_rep1_output_graph.xlsx").dropna()

# T-test between mock and ix for each Target Name

# Take only the data we need
grouped_data = pd.DataFrame(
    {
        "Target Name": data["Target Name"],
        "Sample Name": data["Sample Name"],
        "percent_recovery": data["percent_recovery"],
    }
)
# Empty df to store results of t-test
results = pd.DataFrame(columns=["Target Name", "T-statistic", "P-value"])

# Define target names to use in p-value calcs
target_names = grouped_data["Target Name"].unique()

# Empty list to store p-values for "WT ix" samples
# This way we can graph only these and not superfluous results from all rows
# The number of values will match the number of x axis groups in the final graph
# So we can plot as x = row number
p_values_ix = []

# For each target, get mock and ix values from each row and t-test against each other
# Store results in results df and p-values in the p_values_ix list
for target in target_names:
    target_data = grouped_data[grouped_data["Target Name"] == target]
    mock_values = target_data[target_data["Sample Name"] == "WT mock"][
        "percent_recovery"
    ].values
    ix_values = target_data[target_data["Sample Name"] == "WT ix"][
        "percent_recovery"
    ].values
    t_stat, p_value = ttest_ind(mock_values, ix_values, nan_policy="omit")
    results = results.append(
        {"Target Name": target, "T-statistic": t_stat, "P-value": p_value},
        ignore_index=True,
    )
    # Store WT ix pvalues separately so we can plot them above their columns alone
    if target_data["Sample Name"].str.contains("WT ix").any():
        p_values_ix.append(p_value)
# Merge results with earlier df
grouped_data = grouped_data.merge(results, on="Target Name", how="left").reset_index()

# colors for graph
colors = ["#aeccdb", "#3274a1"]

# Plot mock and ix columns next to each other by hue
g = sns.barplot(
    data=data, x="Target Name", y="percent_recovery", hue="Sample Name", palette=colors
)


# Rotate x-axis labels for all plots
# Adjust the y-axis to 15
# Adjust the margins
plt.gca().set_xticklabels(data["Target Name"].unique(), rotation=45, ha="right")
g.set(ylim=(0, 10))
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.2, top=0.9)
plt.xlabel("")
plt.ylabel("Percent recovery")
plt.title("Cappable-seq on mock and infected cells")

# Add statistical significance annotations for "WT ix" samples
# Define thresholds using if statement to add correct number of stars
# If the p-value is significant add stars where x = Target Name and y = above max recovery value
for i, p_value in enumerate(p_values_ix):
    if p_value < 0.001:
        stars = "***"  # Three stars for p-value < 0.001
    elif p_value < 0.01:
        stars = "**"  # Two stars for p-value < 0.01
    elif p_value < 0.05:
        stars = "*"  # One star for p-value < 0.05
    else:
        stars = ""  # No stars for p-value >= 0.05
    if p_value < 0.05:  # Adjust the significance level as needed
        x = i  # x-coordinate of the bar in p_values_ix
        y = (
            max(data[data["Sample Name"] == "WT ix"]["percent_recovery"]) + 0.5
        )  # y-coordinate above the bar
        plt.annotate(stars, (x + 0.2, y), ha="center", va="bottom", fontsize=12)

# Save the plot
sns.set_style("whitegrid")
plt.savefig("cappable_seq_rep1.svg", format="svg", dpi=300)
