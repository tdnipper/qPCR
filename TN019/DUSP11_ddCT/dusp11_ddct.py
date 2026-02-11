import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# create args to define variables
parser = argparse.ArgumentParser(
    description="Analyze ddCT of a gene when specifying housekeepiag and target genes."
    " Designed to take files from QuantStudio software."
    " Can take replicates if listed under inputs flag."
)

parser.add_argument("-i", "--input", nargs="+", type=str, help="List input qPCR files.")
parser.add_argument("-c", "--control", type=str, help="Control or housekeeping gene.")
parser.add_argument("-t", "--target", type=str, help="Target gene.")
parser.add_argument("-o", "--outfile", type=str, help="Filename for ddCT table.")

args = parser.parse_args()

if args.input:
    infile = args.input

if args.control:
    control = args.control

if args.target:
    target = args.target

# Define functions to remove _ in names and import/format files


def replace_char(s, target_char, replacement_char):
    """Replace a character with another character."""
    return s.replace(target_char, replacement_char)


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
    df["Sample Name"] = df["Sample Name"].apply(
        lambda x: replace_char(x, target_char="_", replacement_char=" ")
    )
    df["Target Name"] = df["Target Name"].apply(
        lambda x: replace_char(x, target_char="_", replacement_char=" ")
    )
    return df


# Use a loop to iterate import_file over all filenames in
files = []
for file in infile:
    file_in = import_file(file)
    files.append(file_in)

concat = pd.concat(files, ignore_index=True)


# Now drop DUSP11KO samples so we only have WT
def drop_row(df: pd.DataFrame, x) -> pd.DataFrame:
    """Drop a row if it contains a string"""
    str_drop = str(x)
    mask = df["Sample Name"].str.contains(str_drop)
    df = df.drop(df[mask].index)
    return df


combo = drop_row(concat, "DUSP11KO")

control_subset = combo.loc[combo["Target Name"].str.contains(control)]
target_subset = combo.loc[combo["Target Name"].str.contains(target)]

# Combine 18s and DUSP11 data, index according to mock or ix, and take delta CT and mean of each
control_subset.index = range(len(control_subset))
control_subset = control_subset.add_suffix(f".{control}")
target_subset.index = range(len(target_subset))
target_subset = target_subset.add_suffix(f".{target}")

data = pd.concat(
    [control_subset, target_subset],
    axis=1,
)
# Use combined data to get deltaCT and clean up dataframe
data["delta"] = data[f"CT.{target}"] - data[f"CT.{control}"]
data = data.drop(
    [
        f"Target Name.{control}",
        f"CT.{control}",
        f"Sample Name.{target}",
        f"Target Name.{target}",
        f"CT.{target}",
    ],
    axis=1,
)
# Use Sample Name column as index to group by mock and ix to get mean and std of deltaCT
data = data.rename(columns={f"Sample Name.{control}": "Sample Name"})
data.set_index("Sample Name", inplace=True)
data.reset_index(inplace=True)

# Separate out mock and ix
data_mock = (
    data.loc[data["Sample Name"].str.contains("mock")].add_suffix(".mock").reset_index()
)

data_ix = (
    data.loc[data["Sample Name"].str.contains("ix")].add_suffix(".ix").reset_index()
)
data_final = pd.concat([data_mock, data_ix], axis=1).drop(
    ["Sample Name.ix", "index"], axis=1
)
# Get control mean (delta.mock.mean) by taking mean of delta.mock
data_final["delta.mock.mean"] = data_final.groupby("Sample Name.mock")[
    "delta.mock"
].transform("mean")
data_final = data_final.rename(
    columns={"Sample Name.mock": "Sample Name"}
).reset_index()
data_final["Sample Name"] = data_final["Sample Name"].apply(
    lambda x: replace_char(x, target_char=" mock", replacement_char="")
)
# Get ddCT for control and target by subtracting delta.ix from delta.mock.mean
data_final[f"ddCT {control}"] = data_final["delta.mock"] - data_final["delta.mock.mean"]
data_final[f"ddCT {target}"] = data_final["delta.ix"] - data_final["delta.mock.mean"]
data_final[f"foldchange {control}"] = np.power(2, -(data_final[f"ddCT {control}"]))
data_final[f"foldchange {target}"] = np.power(2, -(data_final[f"ddCT {target}"]))

melted_df = pd.melt(
    data_final,
    id_vars=["Sample Name"],
    value_vars=["foldchange RNA18S1", "foldchange DUSP11"],
    var_name="foldchange",
)


# Save final table to an excel file
melted_df.to_excel(f"ddCT_{target}_table.xlsx")
