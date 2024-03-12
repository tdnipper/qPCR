#!/usr/bin/env python

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
plt.close()


# Plot melt curve data
def import_melt(filename) -> pd.DataFrame:
    df = pd.read_excel(filename, sheet_name="Melt Curve Raw Data", skiprows=46)
    return df


melt_curve = import_melt("TN019_21.xlsx")

rename = {
    "A1": "mock 16 18S",
    "A2": "mock 16 18S",
    "A3": "mock 16 18S",
    "A4": "mock 16 actin",
    "A5": "mock 16 actin",
    "A6": "mock 16 actin",
    "A7": "mock 16 DUSP11",
    "A8": "mock 16 DUSP11",
    "A9": "mock 16 DUSP11",
    "A10": "infect 16 18S",
    "A11": "infect 16 18S",
    "A12": "infect 16 18S",
    "B1": "infect 16 actin",
    "B2": "infect 16 actin",
    "B3": "infect 16 actin",
    "B4": "infect 16 DUSP11",
    "B5": "infect 16 DUSP11",
    "B6": "infect 16 DUSP11",
    "B7": "mock 24 18S",
    "B8": "mock 24 18S",
    "B9": "mock 24 18S",
    "B10": "mock 24 actin",
    "B11": "mock 24 actin",
    "B12": "mock 24 actin",
    "C1": "mock 24 DUSP11",
    "C2": "mock 24 DUSP11",
    "C3": "mock 24 DUSP11",
    "C4": "infect 24 18S",
    "C5": "infect 24 18S",
    "C6": "infect 24 18S",
    "C7": "infect 24 actin",
    "C8": "infect 24 actin",
    "C9": "infect 24 actin",
    "C10": "infect 24 DUSP11",
    "C11": "infect 24 DUSP11",
    "C12": "infect 24 DUSP11",
}


# Rename columns
def rename_columns(df: pd.DataFrame, rename_dict: dict) -> pd.DataFrame:
    """Rename columns using a dictionary with values for new names."""
    df["Sample Name"] = df["Well Position"].map(rename_dict)
    df = df.dropna()
    return df


melt_curve = rename_columns(melt_curve, rename)
# melt_curve.to_excel("TN019.21 - Melt Curve.xlsx", index=False)


# Plot melt curve
def plot_melt_curve(df: pd.DataFrame, title: str) -> plt:
    df_mock = df[df["Sample Name"].str.contains("mock")]
    df_mock_outlier = df_mock[df_mock["Well Position"].str.contains("A6")]
    df_mock = df_mock[~df_mock["Well Position"].str.contains("A6")]
    df_infect = df[df["Sample Name"].str.contains("infect")]
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 3, 1)
    g_mock = sns.lineplot(
        data=df_mock, x="Temperature", y="Derivative", hue="Sample Name"
    )
    g_mock.set_title(f"Mock - {title}")
    g_mock.set(xlabel="Temperature", ylabel="Derivative")
    plt.subplot(1, 3, 2)
    g_mock_outlier = sns.lineplot(
        data=df_mock_outlier, x="Temperature", y="Derivative", hue="Sample Name"
    )
    g_mock_outlier.set_title(f"Mock Outlier - {title}")
    g_mock_outlier.set(xlabel="Temperature", ylabel="Derivative")
    plt.subplot(1, 3, 3)
    g_infect = sns.lineplot(
        data=df_infect, x="Temperature", y="Derivative", hue="Sample Name"
    )
    g_infect.set_title(f"Infect - {title}")
    g_infect.set(xlabel="Temperature", ylabel="Derivative")
    plt.tight_layout()
    plt.savefig("TN019.21 - Melt Curve.png", format="png", dpi=300)


plot_melt_curve(melt_curve, "Melt Curve")
