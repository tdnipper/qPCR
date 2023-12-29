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


data = import_file("TN037_3.xlsx")


# Sort rows based on input or enrich or flowthrough into new dataframes
def sort_rows(df):
    input_data = pd.DataFrame()
    enrich_data1 = pd.DataFrame()
    enrich_data2 = pd.DataFrame()
    flowthrough1_data = pd.DataFrame()
    flowthrough2_data = pd.DataFrame()
    for index, row in df.iterrows():
        sample_name = row["Sample Name"]
        if "input_beads" in sample_name:
            input_beads_data = pd.concat(
                [input_beads_data, row.to_frame().T], ignore_index=True
            )
        elif "input_no_beads" in sample_name:
            input_no_beads_data = pd.concat(
                [input_no_beads_data, row.to_frame().T], ignore_index=True
            )
        elif "enrich1" in sample_name:
            enrich_data1 = pd.concat(
                [enrich_data1, row.to_frame().T], ignore_index=True
            )
        elif "enrich2" in sample_name:
            enrich_data2 = pd.concat(
                [enrich_data2, row.to_frame().T], ignore_index=True
            )
        elif "flowthrough1_beads" in sample_name:
            flowthrough1_beads_data = pd.concat(
                [flowthrough1_beads_data, row.to_frame().T], ignore_index=True
            )
        elif "flowthrough1_no_beads" in sample_name:
            flowthrough1_no_beads_data = pd.concat(
                [flowthrough1_no_beads_data, row.to_frame().T], ignore_index=True
            )
        elif "flowthrough2" in sample_name:
            flowthrough2_data = pd.concat(
                [flowthrough2_data, row.to_frame().T], ignore_index=True
            )
    return (
        input_beads_data,
        enrich_data1,
        enrich_data2,
        flowthrough1_beads_data,
        flowthrough1_no_beads_data,
        flowthrough2_data,
    )


(
    input_beads_data,
    enrich_data1,
    enrich_data2,
    flowthrough1_beads_data,
    flowthrough1_no_beads_data,
    flowthrough2_data,
) = sort_rows(data)

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

# Add flowthrough data to dataframe for separate graph
flowthrough1 = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
flowthrough1["Sample Name"] = flowthrough1_data["Sample Name"]
flowthrough1["Target Name"] = flowthrough1_data["Target Name"]
flowthrough1["CT.difference"] = flowthrough1_data["CT"].sub(input_data["CT"])
flowthrough1["percent_recovery"] = 100 * 2 ** flowthrough1["CT.difference"]
flowthrough2 = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
flowthrough2["Sample Name"] = flowthrough2_data["Sample Name"]
flowthrough2["Target Name"] = flowthrough2_data["Target Name"]
flowthrough2["CT.difference"] = flowthrough2_data["CT"].sub(input_data["CT"])
flowthrough2["percent_recovery"] = 100 * 2 ** flowthrough2["CT.difference"]

recovery = pd.concat([recover_1, recover_2])
flowthrough = pd.concat([flowthrough1, flowthrough2])

output = pd.concat([recovery, flowthrough])
output.to_excel("TN037.3_output.xlsx")

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
plt.savefig("TN037.3.png", format="png", dpi=300)

f = sns.barplot(
    data=flowthrough,
    x="Target Name",
    y="percent_recovery",
    hue="Sample Name",
    palette="Blues",
    errorbar="sd",
)
