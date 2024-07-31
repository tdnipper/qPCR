import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description= "Calculate ddCT from qPCR data"
)

parser.add_argument(
    "-i", "--infile", type=str, help="File with samples and CTs"
)

parser.add_argument(
    "-o", "--outfile", type=str, help="Name of outfile for ddCT results"
)

args = parser.parse_args()

infile = args.infile
outfile = args.outfile

def import_file(filename, input, control) -> pd.DataFrame:
    """ Import file and parse control and experimental data"""
    df = (
        pd.read_excel(
            filename,
            usecols=["Sample Name", "Target Name", "CT", "Control", "Target"],
            na_values="Undetermined",
        )
    .groupby(by=["Sample Name", "Target Name"])
    )
    control = df["Control"]
    target = df["Target"]
    return df, control, target


    