import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Import function to read results and drop NTC/no amplification
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


data = import_file("TN037_8.xlsx")

# Sort rows based on user defined categories in a given column into new dataframes
def sort_rows(df: pd.DataFrame, categories: list, colname: str) -> tuple[pd.DataFrame]:
    """Sort user defined rows into different user supplied categories and return tuple of dataframes."""
    dfs = {categories: pd.DataFrame() for categories in categories}

    for index, row in df.iterrows():
        sample_name = row[colname]
        for category in categories:
            if category in sample_name:
                dfs[category] = pd.concat(
                    [dfs[category], row.to_frame().T], ignore_index=True
                )
                break
    return tuple(dfs.values())

samples = ["input", "enriched"]

# Get dataframes for each sample to calculate percent recovery
input, enrichment = sort_rows(
    data, samples, "Sample Name"
)

def calculate_percent_recovery(input_df, enrich_df):
    percent_recovery = pd.DataFrame(
        columns=["Sample Name", "Target Name", "CT.difference"]
    )
    percent_recovery["Sample Name"] = enrich_df["Sample Name"]
    percent_recovery["Target Name"] = enrich_df["Target Name"]
    percent_recovery["CT.difference"] = input_df["CT"].sub(enrich_df["CT"])
    percent_recovery["percent_recovery"] = 100 * (2 ** percent_recovery["CT.difference"])
    return percent_recovery

# Calculate percent recovery
enrichment_recovery = calculate_percent_recovery(input, enrichment)
enrichment_recovery["Sample Name"] = enrichment_recovery["Sample Name"].str.split("enrich", expand=True)[0]
 
enrichment_recovery.to_excel("TN037.8.output.xlsx", index=False)

# Plot percent recovery
def plot_percent_recovery(df, title):
    sns.set_theme(style="whitegrid")
    sns.set_palette("deep", 8)
    ax = sns.barplot(x="Sample Name", y="percent_recovery", data=df, hue="Target Name", errorbar="sd")
    ax.set_title(title)
    # ax.set_ylim(0, 100)
    ax.set_xlabel("Target Name")
    ax.set_ylabel("Percent Recovery")
    plt.xticks(rotation=45)
    plt.savefig(f"TN037.8 {title}.png", dpi=300)
    plt.tight_layout(pad=2)
    plt.close()

plot_percent_recovery(enrichment_recovery, "Percent Recovery")