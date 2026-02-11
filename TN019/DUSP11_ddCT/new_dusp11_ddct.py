#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

infile = ["TN019_9.xls", "TN019_15_DUSP11_2.xls", "TN019_18_DUSP11_3.xls"]
files = []
for file in infile:
    file_in = import_file(file)
    files.append(file_in)

concat = pd.concat(files, ignore_index=True)

data = concat[concat["Sample Name"].str.contains("WT")]
t16_mock_18s = data[data["Target Name"].str.contains("18S") & data["Sample Name"].str.contains("16 mock")]
t16_mock_dusp11 = data[data["Target Name"].str.contains("DUSP11") & data["Sample Name"].str.contains("16 mock")]
t16_infect_18s = data[data["Target Name"].str.contains("18S") & data["Sample Name"].str.contains("16 infect")]
t16_infect_dusp11 = data[data["Target Name"].str.contains("DUSP11") & data["Sample Name"].str.contains("16 infect")]
t16_ddCT = pd.DataFrame()
t16_ddCT["Sample Name"] = data["Sample Name"].str.contains("16")
print(t16_ddCT)