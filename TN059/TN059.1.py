#!/usr/bin/env python3

import pandas as pd

def import_file (filename) -> pd.DataFrame:
    """Import the excel file and return a pandas DataFrame."""
    df = pd.read_excel(
            filename,
            sheet_name="Results",
            skiprows=46,
            usecols=["Sample Name", "Target Name", "CT"],
            na_values=["Undetermined", "NTC"]
        ).reset_index(drop=True).dropna()
    return df

data = import_file("TN059_1_redo.xlsx")
data = data[data["Sample Name"].str.contains("noRT") == False].sort_values(by=["Sample Name", "Target Name"]).reset_index(drop=True)

