#!/usr/bin/env python3

import pandas as pd

if __name__ == "__main__":
    print ("This is a class file, import me into a script")

class FileImporter:

    """A class to import qPCR data files formatted specifically from the QuantStudio 7."""
    
    def init (self, filename):
        self.filename = filename

    def import_file (filename) -> pd.DataFrame:
        """Import the excel file and return a pandas DataFrame."""
        df = pd.read_excel(
                filename,
                sheet_name="Results",
                skiprows=46,
                usecols=["Sample Name", "Target Name", "CT"],
                na_values=["Undetermined", "NTC"]
            ).reset_index()
            .dropna()
        )
        return df
    
