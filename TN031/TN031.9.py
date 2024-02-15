import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def import_file(filename):
    """This function will import, take the mean, and format replicate files."""
    df = (
        pd.read_excel(
            filename,
            sheet_name="Results",
            skiprows=46,
            usecols=["Sample Name", "Target Name", "CT"],
            na_values="Undetermined",
        ).dropna()
        # .groupby(by=["Sample Name", "Target Name"])
        # .mean()
        # .reset_index()
    )
    return df


data = import_file("TN031_9.xlsx")


# Sort rows based on user defined categories in a given column into new dataframes
def sort_rows(df, categories, colname):
    """Sort user defined rows into different user supplied categories."""
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


# Sort into 24 and 48 hour timepoints
t24, t48 = sort_rows(data, ["24", "48"], "Sample Name")

# Sort each timepoint into target and housekeeping genes
dusp11_24, actin_24 = sort_rows(t24, ["DUSP11", "actin"], "Target Name")
dusp11_48, actin_48 = sort_rows(t48, ["DUSP11", "actin"], "Target Name")

# Sort each target and housekeeping gene into mock and dox treated
mock_dusp11_24, dox_dusp11_24 = sort_rows(dusp11_24, ["mock", "dox"], "Sample Name")
mock_actin_24, dox_actin_24 = sort_rows(actin_24, ["mock", "dox"], "Sample Name")
mock_dusp11_48, dox_dusp11_48 = sort_rows(dusp11_48, ["mock", "dox"], "Sample Name")
mock_actin_48, dox_actin_48 = sort_rows(actin_48, ["mock", "dox"], "Sample Name")

# Calculate deltaCT between mock and dox treated at each timepoint
dCT_mock = pd.DataFrame()
dCT_dox = pd.DataFrame()
dCT_mock["mock 24"] = mock_dusp11_24["CT"].sub(mock_actin_24["CT"])
dCT_dox["dox 24"] = dox_dusp11_24["CT"].sub(dox_actin_24["CT"])
dCT_mock["mock 48"] = mock_dusp11_48["CT"].sub(mock_actin_48["CT"])
dCT_dox["dox 48"] = dox_dusp11_48["CT"].sub(dox_actin_48["CT"])

ddCT = pd.DataFrame()
ddCT["24 ddCT"] = dCT_dox["dox 24"].sub(dCT_mock["mock 24"])
ddCT["48 ddCT"] = dCT_dox["dox 48"].sub(dCT_mock["mock 48"])

ddCT_foldchange = pd.DataFrame()
ddCT_foldchange["24"] = np.power(2, -(ddCT["24 ddCT"]))
ddCT_foldchange["48"] = np.power(2, -(ddCT["48 ddCT"]))

# Save final table to an excel file
output = pd.concat([ddCT, ddCT_foldchange], axis=1)
output.to_excel("TN031.9 ddCT.xlsx")


# Plot foldchange
def plot_foldchange(df, title):
    """Plot foldchange."""
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 5))
    ax = sns.barplot(data=df, errorbar="se")
    ax.set_title(title)
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("Foldchange")
    plt.savefig(f"TN031.9 {title}.png", dpi=300)


plot_foldchange(ddCT_foldchange, "DUSP11 Foldchange during CRISPRi induction")

# Import melt-curve data

melt_curve = pd.read_excel(
    "TN031_9.xlsx", sheet_name="Melt Curve Raw Data", skiprows=46
)
# Make dictionary to rename columns
rename = {
    "A1": "mock 24 DUSP11",
    "A2": "mock 24 DUSP11",
    "A3": "mock 24 DUSP11",
    "A4": "mock 24 actin",
    "A5": "mock 24 actin",
    "A6": "mock 24 actin",
    "B1": "mock 48 DUSP11",
    "B2": "mock 48 DUSP11",
    "B3": "mock 48 DUSP11",
    "B4": "mock 48 actin",
    "B5": "mock 48 actin",
    "B6": "mock 48 actin",
    "C1": "dox 24 DUSP11",
    "C2": "dox 24 DUSP11",
    "C3": "dox 24 DUSP11",
    "C4": "dox 24 actin",
    "C5": "dox 24 actin",
    "C6": "dox 24 actin",
    "D1": "dox 48 DUSP11",
    "D2": "dox 48 DUSP11",
    "D3": "dox 48 DUSP11",
    "D4": "dox 48 actin",
    "D5": "dox 48 actin",
    "D6": "dox 48 actin",
}


# Rename columns
def rename_columns(df, rename):
    """Rename columns."""
    df["Sample Name"] = df["Well Position"].map(rename)
    df = df.dropna()
    return df


melt_curve = rename_columns(melt_curve, rename)


# Plot melt curve
def plot_melt_curve(df, title):
    """Plot melt curve."""
    # Separate into mock and dox treated dataframes
    df_mock = df[df["Sample Name"].str.contains("mock")]
    df_dox = df[df["Sample Name"].str.contains("dox")]
    # Combine mock and dox plots into one figure of separate plots
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    g_mock = sns.lineplot(
        data=df_mock, x="Temperature", y="Derivative", hue="Sample Name"
    )
    g_mock.set_title(f"Mock - {title}")
    g_mock.set(xlabel="Temperature", ylabel="Derivative")
    plt.subplot(1, 2, 2)
    g_dox = sns.lineplot(
        data=df_dox, x="Temperature", y="Derivative", hue="Sample Name"
    )
    g_dox.set_title(f"Dox Treated - {title}")
    g_dox.set(xlabel="Temperature", ylabel="Derivative")

    plt.tight_layout()
    plt.savefig(f"TN031.9 Combined - {title}.png", dpi=300)


plot_melt_curve(melt_curve, "Melt Curve")
