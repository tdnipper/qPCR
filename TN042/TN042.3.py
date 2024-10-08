#!/usr/bin/env python3

import pandas as pd
import altair as alt
import plotly.express as px

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


data = import_file("TN042_3.xlsx")


# Filter the data based on target gene and take the mean of the cT for each Sample Name
def filter_targets(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter the data based on the target gene and take the mean of the cT."""
    filtered = df[df["Target Name"].str.contains(target)].copy()
    filtered["mean"] = filtered.groupby("Sample Name")["CT"].transform("mean")
    return filtered


target_dict = {}
targets = ["18S", "DUSP11", "IFNB", "MX1"] # These must match the input data including case
for target in targets:
    target_dict[target.lower()] = filter_targets(data, target)


# Sort data by sample name (+- dox and ix) and remove duplicates, keeping only CT mean
def process_data(df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    return (
        df[df["Sample Name"] == sample_name]
        .drop("CT", axis=1)
        .drop_duplicates(subset=["Sample Name", "mean"])
    )


sample_names = ["+dox mock", "-dox mock", "+dox ix", "-dox ix"]
results = {}

for key, data in target_dict.items():
    for sample in sample_names:
        results[f"{key}_{sample}"] = process_data(data, sample)


# Calculate ddCT during infection for all targets during both dox treatment states
def calculate_ddCT_ix(results, gene, treat, control_gene="18s"):
    return (
        results[f"{gene}_{treat} ix"]["mean"].values
        - results[f"{control_gene}_{treat} ix"]["mean"].values
    ) - (
        results[f"{gene}_{treat} mock"]["mean"].values
        - results[f"{control_gene}_{treat} mock"]["mean"].values
    )


# Calculate ddCT during dox treatment for mock samples only
def calculate_ddCT_dox(results, gene, control_gene="18s"):
    return (
        results[f"{gene}_+dox mock"]["mean"].values
        - results[f"{control_gene}_+dox mock"]["mean"].values
    ) - (
        results[f"{gene}_-dox mock"]["mean"].values
        - results[f"{control_gene}_-dox mock"]["mean"].values
    )

# Note: can probably combine these two functions if we define the dox and ix states

ddct_results = {}
ddct_results_dox = {}
genes = ["dusp11", "ifnb", "mx1"]
treatments = ["-dox", "+dox"]

for treat in treatments:
    for gene in genes:
        ddct_results[f"{gene}_{treat}"] = calculate_ddCT_ix(results, gene, treat)

for gene in genes:
    ddct_results_dox[gene] = calculate_ddCT_dox(results, gene)

# Generate ddCT and foldchange dataframes during infection for all conditions
ddct_df_ix = pd.DataFrame()
foldchange_df_ix = pd.DataFrame()
for key, ddct in ddct_results.items():
    ddct_df_ix[key] = ddct
    foldchange_df_ix[f"{key} foldchange"] = 2**-ddct
foldchange_df_ix = foldchange_df_ix.melt(var_name="Condition", value_name="Fold Change")
ddct_df_ix = ddct_df_ix.melt(var_name="Condition", value_name="ddCT")
foldchange_df_ix["dox"] = foldchange_df_ix["Condition"].str.split("_").str[1]
foldchange_df_ix["dox"] = foldchange_df_ix["dox"].str.split(" ").str[0]
foldchange_df_ix["Condition"] = foldchange_df_ix["Condition"].str.split("_").str[0]
ddct_df_ix.to_excel("ddct_ix.xlsx")
foldchange_df_ix.to_excel("foldchange_ix.xlsx")

# Generate ddCT and foldchange dataframes during dox induction only for mock infected cells
ddct_df_dox = pd.DataFrame()
foldchange_df_dox = pd.DataFrame()
for key, ddct in ddct_results_dox.items():
    ddct_df_dox[key] = ddct
    foldchange_df_dox[f"{key} foldchange"] = 2**-ddct
ddct_df_dox = ddct_df_dox.melt(var_name="Condition", value_name="ddCT")
ddct_df_dox.to_excel("ddct_dox.xlsx")
foldchange_df_dox = foldchange_df_dox.melt(var_name="Condition", value_name="Fold Change")
foldchange_df_dox["Condition"] = foldchange_df_dox["Condition"].str.split("foldchange").str[0]
foldchange_df_dox.to_excel("foldchange_dox.xlsx")

# Try altair out
title = alt.TitleParams("Fold Change During Infection", anchor="middle")
ix_plot = alt.Chart(foldchange_df_ix, title = title).mark_bar().encode(
    alt.X("Condition", title="", axis=alt.Axis(labelAngle=0)),
    alt.Y("Fold Change", title="Fold Change vs. Mock Infected"),
    xOffset=alt.XOffset("dox:N"),
    color = alt.Color("dox").title("Dox Treatment"),
).properties(width=300)
ix_plot.save("ix_plot.png", ppi=300)

dox_plot = alt.Chart(foldchange_df_dox, title = "Fold Change During Dox Treatment").mark_bar().encode(
    alt.X("Condition", title="", axis=alt.Axis(labelAngle=0)),
    alt.Y("Fold Change", title="Fold Change vs. Untreated")
).properties(width=300)
dox_plot.save("dox_plot.png", ppi=300)

# Try plotly out, much shorter than altair and graphs are interactive when html
fig = px.bar(foldchange_df_ix, x="Condition", y="Fold Change", color="dox", barmode="group", title="Fold Change During Infection")
fig.write_image("ix_plot_plotly.svg")