import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

cwd=os.getcwd()
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


data = import_file("TN037_2.xlsx")


# Sort rows based on input or enrich into new dataframes
def sort_rows(df):
    input_data = pd.DataFrame()
    enrich_data1 = pd.DataFrame()
    enrich_data2 = pd.DataFrame()
    for index, row in df.iterrows():
        sample_name = row["Sample Name"]
        if "input" in sample_name:
            input_data = pd.concat([input_data, row.to_frame().T], ignore_index=True)
        elif "enriched 1" in sample_name:
            enrich_data1 = pd.concat(
                [enrich_data1, row.to_frame().T], ignore_index=True
            )
        elif "enriched 2" in sample_name:
            enrich_data2 = pd.concat(
                [enrich_data2, row.to_frame().T], ignore_index=True
            )
    return input_data, enrich_data1, enrich_data2


input_data, enrich_data1, enrich_data2 = sort_rows(data)

# Calculate difference between input and each enrichment CTs
# Use this number to calculate percent recovery
recover_1 = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover_1["Sample Name"] = enrich_data1["Sample Name"]
recover_1["Target Name"] = enrich_data1["Target Name"]
recover_1["CT.difference"] = input_data["CT"].sub(enrich_data1["CT"])
recover_1["percent_recovery"] = 100 * 2 ** recover_1["CT.difference"]
recover_2 = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover_2["Sample Name"] = enrich_data2["Sample Name"]
recover_2["Target Name"] = enrich_data2["Target Name"]
recover_2["CT.difference"] = input_data["CT"].sub(enrich_data2["CT"])
recover_2["percent_recovery"] = 100 * 2 ** recover_2["CT.difference"]
recovery = pd.concat([recover_1, recover_2])

recovery.to_excel("TN037.2_output.xlsx")

# Plot using paried barplot to compare both enrichments
g = sns.barplot(
    data=recovery,
    x="Target Name",
    y="percent_recovery",
    hue="Sample Name",
    palette="Blues",
    errorbar="sd",
)
plt.legend(title="Enrichment")
plt.ylabel("Percent recovery relative to unenriched")
plt.savefig("TN037.2.png", format="png", dpi=300)
