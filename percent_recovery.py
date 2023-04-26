import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Analyse % recovery from cappable seq files"
)
parser.add_argument(
    "-i", "--input", nargs="+", type=str, help="File with input samples"
)
parser.add_argument(
    "-k",
    "--key",
    type=str,
    help="Key file with names of transcripts to be analyzed in a column",
)

parser.add_argument("-e", "--enrich", type=str, help="File with enriched samples")

args = parser.parse_args()
infile = args.input
efile = args.enrich
if args.key:
    key = args.key

# Get the keyfile provided and specify the column that names are in
# keyfile = pd.read_excel(key)
# keycol_name = keyfile.columns[1]

# # Create a list of genes as a dictionary from the keyfile provided by the user
# gene_dict = {}
# for i, val in enumerate(keyfile[keycol_name]):
#     gene_dict[f"var{i}"] = val


# Import files for subsequent analysis
def import_file(filename) -> pd.DataFrame:
    df = pd.read_excel(filename)


# Use a loop to iterate import_file over all filenames in
files_input = []
for file in infile:
    file_in = import_file(file)
    files_input.append(file_in)
files_enrich = []
for file in efile:
    file_in = import_file(file)
    files_enrich.append(file_in)

concat_input = pd.concat(files_input, ignore_index=True)
concat_enrich = pd.concat(files_enrich, ignore_index=True)

print(concat_enrich)
