import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


data = import_file("TN043_5.xlsx")

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
foldchange_df_ix_dusp11 = foldchange_df_ix[foldchange_df_ix["Sample"].str.contains("dusp11")].copy()
# Normalize DUSP11 fold change to T0
foldchange_df_ix_dusp11.loc[:,"normalized"] = foldchange_df_ix_dusp11["Fold Change"]/foldchange_df_ix_dusp11["Fold Change"].max()
foldchange_df_ix_ifnb = foldchange_df_ix[foldchange_df_ix["Sample"].str.contains("ifnb")]

# Plot the data
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
fig1 = sns.barplot(
    data = foldchange_df_ix_dusp11,
    x = "Sample",
    y = "normalized",
    ax=ax[0],

)

fig2 = sns.barplot(
    data = foldchange_df_ix_ifnb,
    x = "Sample",
    y = "Fold Change",
    ax=ax[1]
)
ax[0].set_xticks(range(6))
ax[1].set_xticks(range(6))
labels = ["0", "1", "2", "4", "6", "8"]
ax[0].set_xticklabels(labels)
ax[1].set_xticklabels(labels)
ax[0].set_title("DUSP11 mRNA Fold Change")
ax[1].set_title("IFNB mRNA Fold Change")
ax[0].set_ylabel("Fold Change")
ax[1].set_ylabel("")
ax[0].set_xlabel("Hours post-infection")
ax[1].set_xlabel("Hours post-infection")
plt.tight_layout()
plt.savefig("TN043_5.png", dpi=300)