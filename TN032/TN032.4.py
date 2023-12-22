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
        .dropna()
    )
    return df


data = import_file("TN032_4_redo.xlsx")

# Extract relevant subsets
genes_of_interest = ["RNA18S1", "WSN_NP", "DUSP11"]  # used to filter data
subset_data = (
    data[data["Target Name"].str.contains("|".join(genes_of_interest))].pivot(
        index="Sample Name", columns="Target Name", values="CT"
    )  # re organize so that Target Names are now separated into columns using the sample name as the index
    # values will be CT
    .reset_index()
)

# Calculate delta values
subset_data["delta.dusp11"] = subset_data["DUSP11"] - subset_data["RNA18S1"]
subset_data["delta.NP"] = subset_data["WSN_NP"] - subset_data["RNA18S1"]
# Separate into mock and infected samples
mock = subset_data[subset_data["Sample Name"].str.contains("T0")].reset_index(drop=True)
mock_dusp11 = mock.at[0, "delta.dusp11"]
mock_NP = 40
ix = subset_data[subset_data["Sample Name"].str.contains("T")].reset_index(drop=True)
ix["Sample Name"] = ix["Sample Name"].str.extract("(\d+)").astype(int)
ix = ix.sort_values(by="Sample Name").reset_index(drop=True)

# Calculate ddCT
ddct = pd.DataFrame()
ddct["Name"] = ["T0", "T1", "T2", "T4", "T7", "T24"]
ddct["dusp11"] = ix["delta.dusp11"].subtract(mock_dusp11)
ddct["NP"] = ix["delta.NP"].subtract(mock_NP)
ddct_melted = pd.melt(ddct, id_vars=["Name"], var_name="Gene", value_name="CT Change")
ddct_melted["foldchange"] = np.power(2, -(ddct_melted["CT Change"]))

plot = sns.barplot(
    data=ddct_melted, x="Name", y="foldchange", hue="Gene", palette="Paired"
)
plt.title("DUSP11 mRNA fold change during infection")
plt.ylabel("fold change")
plt.savefig("TN032.4_redo.png", format="png", dpi=300)
plt.ylim(0, 1.2)
plt.savefig("TN032.4_redo_zoomed.png", format="png", dpi=300)
ddct_melted.to_excel("TN032.4_redo_output.xlsx")
