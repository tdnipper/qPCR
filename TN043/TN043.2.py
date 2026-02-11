#!/usr/bin/env python3

import pandas as pd
# import altair as alt
# import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt

def import_file(filename) -> pd.DataFrame:
    """This function will import, take the mean, and format replicate files."""
    df = (
        pd.read_excel(
            filename,
            sheet_name="Results",
            skiprows=46,
            usecols=["Sample Name", "Target Name", "CT"],
            na_values=["Undetermined", "NTC"],
        )
        .reset_index()
        .dropna()
    )
    return df


data = import_file("TN043_2.xlsx")

# Filter the data based on target gene and take the mean of the cT for each Sample Name
def filter_targets(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter the data based on the target gene and take the mean of the cT."""
    filtered = df[df["Target Name"].str.contains(target)].copy()
    filtered["mean"] = filtered.groupby("Sample Name")["CT"].transform("mean")
    return filtered

# Sort data by sample name (+- dox and ix) and remove duplicates, keeping only CT mean
def process_data(df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    return (
        df[df["Sample Name"] == sample_name]
        .drop("CT", axis=1)
        .drop_duplicates(subset=["Sample Name", "mean"])
    )


target_dict = {}
targets = ["RNA18S1", "DUSP11", "IFNB"] # These must match the input data including case
for target in targets:
    target_dict[target.lower()] = filter_targets(data, target)

sample_names = ["mock", "T0", "T1", "T2", "T4", "T6", "T8"]
results = {}

for key, data in target_dict.items():
    for sample in sample_names:
        results[f"{key}_{sample}"] = process_data(data, sample)

# Calculate ddCT during infection for all targets during infection relative to mock
def calculate_ddCT(results, gene, timepoint, control_gene="rna18s1", control_timepoint="mock"):
    return (
        results[f"{gene}_{timepoint}"]["mean"].values
        - results[f"{control_gene}_{timepoint}"]["mean"].values
    ) - (
        results[f"{gene}_{control_timepoint}"]["mean"].values
        - results[f"{control_gene}_{control_timepoint}"]["mean"].values
    )

ddct_results = {}
genes = ["dusp11", "ifnb"]
timepoints = ["T0", "T1", "T2", "T4", "T6", "T8"]

for gene in genes:
    for timepoint in timepoints:
        ddct_results[f"{gene}_{timepoint}"] = calculate_ddCT(results, gene, timepoint)

# Calculate foldchange for all targets during infection relative to mock
ddct_df_ix = pd.DataFrame()
foldchange_df_ix = pd.DataFrame()
for key, ddct in ddct_results.items():
    ddct_df_ix[key] = ddct
    foldchange_df_ix[f"{key} foldchange"] = 2**-ddct
foldchange_df_ix = foldchange_df_ix.melt(var_name="Sample", value_name="Fold Change")

# # Plot the fold change during infection
# fig = px.bar(foldchange_df_ix[foldchange_df_ix["Sample"].str.contains("dusp11")], x="Sample", y="Fold Change", title="Fold Change During Infection")
# fig.write_image("TN043.2_dusp11.svg")
# fig.write_image("TN043.2_dusp11.png")
# fig2 = px.bar(foldchange_df_ix[foldchange_df_ix["Sample"].str.contains("ifnb")], x="Sample", y="Fold Change", title="Fold Change During Infection")
# fig2.write_image("TN043.2_ifnb.svg")
# fig2.write_image("TN043.2_ifnb.png")

fig = sns.barplot(data=foldchange_df_ix[foldchange_df_ix["Sample"].str.contains("dusp11")], x="Sample", y="Fold Change")
plt.title("DUSP11 Fold Change During Infection")
labels = ["T0", "T1", "T2", "T4", "T6", "T8"]
plt.xticks(range(6), labels)
plt.xlabel("")
plt.tight_layout()
plt.savefig("TN043.2_dusp11.png", dpi=300)