import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Analyse % recovery from cappable seq files"
)
parser.add_argument(
    "-i", "--input", nargs="+", type=str, help="List of input excel files"
)
parser.add_argument(
    "-k",
    "--key",
    type=str,
    help="Key file with names of transcripts to be analyzed in a column",
)

args = parser.parse_args()
infile = args.input
key = args.key

# def import_key(filename) -> pd.DataFrame:
#     pd.read_excel(filename)

# Get the keyfile provided and specify the column that names are in
keyfile = pd.read_excel(key)
keycol_name = keyfile.columns[1]

# Create a list of genes as a dictionary from the keyfile provided by the user
gene_dict = {}
for i, val in enumerate(keyfile[keycol_name]):
    gene_dict[f"var{i}"] = val


# Import files for subsequent analysis
def import_file(filename) -> pd.DataFrame:
    df = pd.read_excel(filename)


# Use a loop to iterate import_file over all filenames in
files = []
for file in infile:
    file_in = import_file(file)
    files.append(file_in)

concat = pd.concat(files, ignore_index=True)
