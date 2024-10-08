import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
os.chdir(cwd)


# Import function to read results and take the mean of technical replicates
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


data = import_file("TN037_5.xlsx")


# Sort rows based on input or enrich or flowthrough into new dataframes
def sort_rows(df, categories):
    dfs = {categories: pd.DataFrame() for categories in categories}

    for index, row in df.iterrows():
        sample_name = row["Sample Name"]
        for category in categories:
            if category in sample_name:
                dfs[category] = pd.concat(
                    [dfs[category], row.to_frame().T], ignore_index=True
                )
                break
    return tuple(dfs.values())

samples = list(data["Sample Name"].unique())

# Get dataframes for each sample to calculate percent recovery
Zymo_input, Trizol_input, Zymo_old, Trizol_old, Zymo_new, Trizol_new = sort_rows(
    data, samples
)


def calculate_percent_recovery(input_df, enrich_df):
    percent_recovery = pd.DataFrame(
        columns=["Sample Name", "Target Name", "CT.difference"]
    )
    percent_recovery["Sample Name"] = enrich_df["Sample Name"]
    percent_recovery["Target Name"] = enrich_df["Target Name"]
    percent_recovery["CT.difference"] = input_df["CT"].sub(enrich_df["CT"])
    percent_recovery["percent_recovery"] = 100 * 2 ** percent_recovery["CT.difference"]
    return percent_recovery


# Calculate percent recovery
recover_zymo_old = calculate_percent_recovery(Zymo_input, Zymo_old)

recover_trizol_old = calculate_percent_recovery(Trizol_input, Trizol_old)

recover_zymo_new = calculate_percent_recovery(Zymo_input, Zymo_new)

recover_trizol_new = calculate_percent_recovery(Trizol_input, Trizol_new)

input_compare = calculate_percent_recovery(Trizol_input, Zymo_input)

zymo = pd.concat([recover_zymo_old, recover_zymo_new])
trizol = pd.concat([recover_trizol_old, recover_trizol_new])
output = pd.concat([zymo, trizol])
output.to_excel("TN037.5_output.xlsx")

# Plot percent recovery
def plot_percent_recovery(df, title):
    sns.set_theme(style="whitegrid")
    sns.set_palette("deep", 8)
    ax = sns.barplot(x="Target Name", y="percent_recovery", data=df, hue="Sample Name")
    ax.set_title(title)
    # ax.set_ylim(0, 100)
    ax.set_xlabel("Target Name")
    ax.set_ylabel("Percent Recovery")
    plt.savefig(f"TN037.5 {title}.png", dpi=300)
    plt.close()


plot_percent_recovery(zymo, "Zymo Recovery")
plot_percent_recovery(trizol, "Trizol Recovery")
plot_percent_recovery(input_compare, "Input Comparison Trizol vs Zymo")
print('Done!')
