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


data = import_file("TN037_4.xlsx")


# Sort rows based on input or enrich or flowthrough into new dataframes
def sort_rows(df):
    input_old = pd.DataFrame()
    input_new = pd.DataFrame()
    flowthrough1_old = pd.DataFrame()
    flowthrough1_new = pd.DataFrame()
    enrich1_old = pd.DataFrame()
    enrich1_new = pd.DataFrame()
    enrich2_old = pd.DataFrame()
    enrich2_new = pd.DataFrame()
    for index, row in df.iterrows():
        sample_name = row["Sample Name"]
        if "input old" in sample_name:
            input_old = pd.concat([input_old, row.to_frame().T], ignore_index=True)
        elif "input new" in sample_name:
            input_new = pd.concat([input_new, row.to_frame().T], ignore_index=True)
        elif "FT old" in sample_name:
            flowthrough1_old = pd.concat(
                [flowthrough1_old, row.to_frame().T], ignore_index=True
            )
        elif "FT new" in sample_name:
            flowthrough1_new = pd.concat(
                [flowthrough1_new, row.to_frame().T], ignore_index=True
            )
        elif "E1 old" in sample_name:
            enrich1_old = pd.concat([enrich1_old, row.to_frame().T], ignore_index=True)
        elif "E1 new" in sample_name:
            enrich1_new = pd.concat([enrich1_new, row.to_frame().T], ignore_index=True)
        elif "E2 old" in sample_name:
            enrich2_old = pd.concat([enrich2_old, row.to_frame().T], ignore_index=True)
        elif "E2 new" in sample_name:
            enrich2_new = pd.concat([enrich2_new, row.to_frame().T], ignore_index=True)
    return (
        input_old,
        input_new,
        flowthrough1_old,
        flowthrough1_new,
        enrich1_old,
        enrich1_new,
        enrich2_old,
        enrich2_new,
    )


(
    input_old,
    input_new,
    flowthrough1_old,
    flowthrough1_new,
    enrich1_old,
    enrich1_new,
    enrich2_old,
    enrich2_new,
) = sort_rows(data)

# Calculate percent recovery for each sample

recover1_old = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover1_old["Sample Name"] = enrich1_old["Sample Name"]
recover1_old["Target Name"] = enrich1_old["Target Name"]
recover1_old["CT.difference"] = input_old["CT"].sub(enrich1_old["CT"])
recover1_old["percent_recovery"] = 100 * 2 ** recover1_old["CT.difference"]

recover1_new = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover1_new["Sample Name"] = enrich1_new["Sample Name"]
recover1_new["Target Name"] = enrich1_new["Target Name"]
recover1_new["CT.difference"] = input_new["CT"].sub(enrich1_new["CT"])
recover1_new["percent_recovery"] = 100 * 2 ** recover1_new["CT.difference"]

recover2_old = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover2_old["Sample Name"] = enrich2_old["Sample Name"]
recover2_old["Target Name"] = enrich2_old["Target Name"]
recover2_old["CT.difference"] = input_old["CT"].sub(enrich2_old["CT"])
recover2_old["percent_recovery"] = 100 * 2 ** recover2_old["CT.difference"]

recover2_new = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
recover2_new["Sample Name"] = enrich2_new["Sample Name"]
recover2_new["Target Name"] = enrich2_new["Target Name"]
recover2_new["CT.difference"] = input_new["CT"].sub(enrich2_new["CT"])
recover2_new["percent_recovery"] = 100 * 2 ** recover2_new["CT.difference"]

flowthrough_old = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
flowthrough_old["Sample Name"] = flowthrough1_old["Sample Name"]
flowthrough_old["Target Name"] = flowthrough1_old["Target Name"]
flowthrough_old["CT.difference"] = input_old["CT"].sub(flowthrough1_old["CT"])
flowthrough_old["percent_recovery"] = 100 * 2 ** flowthrough_old["CT.difference"]
flowthrough_old["percent_recovery_corrected"] = (
    100 * 2 ** flowthrough_old["CT.difference"]
) / 5

flowthrough_new = pd.DataFrame(columns=["Sample Name", "Target Name", "CT.difference"])
flowthrough_new["Sample Name"] = flowthrough1_new["Sample Name"]
flowthrough_new["Target Name"] = flowthrough1_new["Target Name"]
flowthrough_new["CT.difference"] = input_new["CT"].sub(flowthrough1_new["CT"])
flowthrough_new["percent_recovery"] = 100 * 2 ** flowthrough_new["CT.difference"]
flowthrough_new["percent_recovery_corrected"] = (
    100 * 2 ** flowthrough_new["CT.difference"]
) / 5

# recover_old = pd.concat([recover1_old, recover2_old])
# recover_new = pd.concat([recover1_new, recover2_new])
recover = pd.concat([recover1_old, recover2_old, recover1_new, recover2_new])
flowthrough = pd.concat([flowthrough_old, flowthrough_new])

output = pd.concat([recover_old, recover_new, flowthrough])
output.to_excel("TN037.4_output.xlsx")

# Plot using paried barplot to compare both enrichments
g = sns.barplot(
    data=recover_old,
    x="Target Name",
    y="percent_recovery",
    hue="Sample Name",
    palette="Blues",
    errorbar="sd",
)
g.set_title("Enrichment")
g.set_ylabel("Percent recovery relative to input")
plt.tight_layout()
g.figure.savefig("TN037.4_old.png", format="png", dpi=300)
plt.close()

# f = sns.barplot(
#     data=recover_new,
#     x="Target Name",
#     y="percent_recovery",
#     hue="Sample Name",
#     palette="Blues",
#     errorbar="sd",
# )
# f.set_title("Enrichment")
# f.set_ylabel("Percent recovery relative to input")
# plt.tight_layout()
# f.figure.savefig("TN037.4_new.png", format="png", dpi=300)
# plt.close()

# h = sns.barplot(
#     data=flowthrough,
#     x="Target Name",
#     y="percent_recovery",
#     hue="Sample Name",
#     palette="Blues",
#     errorbar="sd",
# )
# h.set_title("Flowthrough")
# h.set_ylabel("Percent recovery relative to input")
# plt.tight_layout()
# h.figure.savefig("TN037.4_flowthrough.png", format="png", dpi=300)
# plt.close()

i = sns.barplot(
    data=flowthrough,
    x="Target Name",
    y="percent_recovery_corrected",
    hue="Sample Name",
    palette="Blues",
    errorbar="sd",
)
i.set_title("Flowthrough Corrected")
i.set_ylabel("Percent recovery relative to input")
plt.tight_layout()
i.figure.savefig("TN037.4_flowthrough_corrected.png", format="png", dpi=300)
plt.close()

j = sns.barplot(
    data=recover,
    x="Target Name",
    y="percent_recovery",
    hue="Sample Name",
    palette="Paired",
    errorbar="sd",
)
j.set_title("Enrichment")
j.set_ylabel("Percent recovery relative to input")
plt.tight_layout()
j.figure.savefig("TN037.4_combined.png", format="png", dpi=300)
plt.close()
