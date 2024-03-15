import pandas as pd


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


# Create separate dataframes using 18S and target genes
def filter_targets(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter the data based on the target gene and take the mean of the cT."""
    filtered = df[df["Target Name"].str.contains(target)].copy()
    filtered["mean"] = filtered.groupby("Sample Name")["CT"].transform("mean")
    return filtered


# Filter data by target gene
target_dict = {}
targets = ["18S", "DUSP11", "IFNB", "MX1"]
for target in targets:
    target_dict[target] = filter_targets(data, target)

# Define function to sort data by sample name and remove duplicates
def process_data(df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    return (
        df[df["Sample Name"] == sample_name]
        .drop("CT", axis=1)
        .drop_duplicates(subset=["Sample Name", "mean"])
    )


data_dict = {
    "18S": target_dict["18S"],
    "dusp11": target_dict["DUSP11"],
    "ifnb": target_dict["IFNB"],
    "mx1": target_dict["MX1"],
}
sample_names = ["+dox mock", "-dox mock", "+dox ix", "-dox ix"]
results = {}

# Process data for each target gene and sample name
for key, data in data_dict.items():
    for sample in sample_names:
        results[f"{key}_{sample}"] = process_data(data, sample)


# Calculate ddCT during infection for all targets during both dox treatment states
def calculate_ddCT_ix(results, gene, treat, control_gene="18S"):
    return (
        results[f"{gene}_{treat} ix"]["mean"].values
        - results[f"{control_gene}_{treat} ix"]["mean"].values
    ) - (
        results[f"{gene}_{treat} mock"]["mean"].values
        - results[f"{control_gene}_{treat} mock"]["mean"].values
    )


# Calculate ddCT during dox treatment for mock samples only
def calculate_ddCT_dox(results, gene, control_gene="18S"):
    return (
        results[f"{gene}_+dox mock"]["mean"].values
        - results[f"{control_gene}_+dox mock"]["mean"].values
    ) - (
        results[f"{gene}_-dox mock"]["mean"].values
        - results[f"{control_gene}_-dox mock"]["mean"].values
    )


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
ddct_df_ix.to_excel("ddct_ix.xlsx")
foldchange_df_ix.to_excel("foldchange_ix.xlsx")

# Generate ddCT and foldchange dataframes during dox induction only for mock infected cells
ddct_df_dox = pd.DataFrame()
foldchange_df_dox = pd.DataFrame()
for key, ddct in ddct_results_dox.items():
    ddct_df_dox[key] = ddct
    foldchange_df_dox[f"{key} foldchange"] = 2**-ddct
ddct_df_dox.to_excel("ddct_dox.xlsx")
foldchange_df_dox.to_excel("foldchange_dox.xlsx")
