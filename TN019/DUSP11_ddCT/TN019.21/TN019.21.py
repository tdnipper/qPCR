import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Import function to read results and take the mean of technical replicates
def import_file(filename) -> pd.DataFrame:
    """This function will import, take the mean, and format replicate files."""
    df = (
        pd.read_excel(
            filename,
            sheet_name="Results",
            skiprows=46,
            usecols=["Sample Name", "Target Name", "CT"],
            na_values="Undetermined",
        )
        .groupby(by=["Sample Name", "Target Name"])
        .mean()
        .reset_index()
    )
    return df


data = import_file("TN019_21.xlsx")

# Extract relevant subsets
genes_of_interest = ["RNA18S1", "B-actin", "DUSP11"]  # used to filter data
subset_data = (
    data[data["Target Name"].str.contains("|".join(genes_of_interest))].pivot(
        index="Sample Name", columns="Target Name", values="CT"
    )  # re organize so that Target Names are now separated into columns using the sample name as the index
    # values will be CT
    .reset_index()
)


# Calculate delta values
subset_data["delta.dusp11"] = subset_data["DUSP11"] - subset_data["RNA18S1"]
subset_data["delta.actin"] = subset_data["B-actin"] - subset_data["RNA18S1"]

# Separate into mock and infected samples
mock = subset_data[subset_data["Sample Name"].str.contains("mock")].reset_index(
    drop=True
)
ix = subset_data[subset_data["Sample Name"].str.contains("ix")].reset_index(drop=True)

# Calculate ddCT
ddct = pd.DataFrame()
ddct["Name"] = ["T16", "T24"]
ddct["dusp11"] = ix["delta.dusp11"] - mock["delta.dusp11"]
ddct["actin"] = ix["delta.actin"] - mock["delta.actin"]
ddct_melted = pd.melt(ddct, id_vars=["Name"], var_name="Gene", value_name="CT Value")
ddct_melted["foldchange"] = np.power(2, -(ddct_melted["CT Value"]))

plot = sns.barplot(
    data=ddct_melted, x="Name", y="foldchange", hue="Gene", palette="Paired"
)
plt.savefig("TN019.21.png", format="png", dpi=300)
plt.savefig("TN019.21.svg", format="svg", dpi=300)