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


# Create control dataframes using 18S and target genes
def filter_targets(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter the data based on the target gene and take the mean of the cT."""
    filtered = df[df["Target Name"].str.contains(target)].copy()
    filtered["mean"] = filtered.groupby("Sample Name")["CT"].transform("mean")
    return filtered


data_18S = filter_targets(data, "18S")
data_dusp11 = filter_targets(data, "DUSP11")
data_ifnb = filter_targets(data, "IFNB")
data_mx1 = filter_targets(data, "MX1")

# Make mock sorted dataframes


def process_data(df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    return (
        df[df["Sample Name"] == sample_name]
        .drop("CT", axis=1)
        .drop_duplicates(subset=["Sample Name", "mean"])
    )


data_dict = {"18S": data_18S, "dusp11": data_dusp11, "ifnb": data_ifnb, "mx1": data_mx1}

sample_names = ["+dox mock", "-dox mock", "+dox ix", "-dox ix"]

results = {}

for key, data in data_dict.items():
    for sample in sample_names:
        results[f"{key}_{sample}"] = process_data(data, sample)
dCT_dusp11_dox = (
    results["dusp11_+dox mock"]["mean"].values - results["18S_+dox mock"]["mean"].values
)

# Calculate ddCT during infection for all targets during both dox treatment states
def calculate_ddCT(results, gene, treat, control_gene="18S"):
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
        ddct_results[f"{gene}_{treat}"] = calculate_ddCT(results, gene, treat)

for gene in genes:
    ddct_results_dox[gene] = calculate_ddCT_dox(results, gene)

ddct_df_ix = pd.DataFrame()
foldchange_df_ix = pd.DataFrame()
for key, ddct in ddct_results.items():
    ddct_df_ix[key] = ddct
    foldchange_df_ix[f'{key} foldchange'] = 2 ** -ddct
ddct_df_ix.to_excel("ddct_ix.xlsx")
foldchange_df_ix.to_excel("foldchange_ix.xlsx")

ddct_df_dox = pd.DataFrame()
foldchange_df_dox = pd.DataFrame()    
for key, ddct in ddct_results_dox.items():
    ddct_df_dox[key] = ddct
    foldchange_df_dox[f"{key} foldchange"] = 2 ** -ddct
ddct_df_dox.to_excel("ddct_dox.xlsx")
foldchange_df_dox.to_excel("foldchange_dox.xlsx")
