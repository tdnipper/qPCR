import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="Analyse % recovery from cappable seq files"
)
parser.add_argument(
    "-i", "--input", nargs="+", type=str, help="File with input samples"
)

parser.add_argument(
    "-e", "--enrich", type=str, nargs="+", help="File with enriched samples"
)

parser.add_argument(
    "-o", "--outfile", type=str, help="Name of outfile for graph and table results"
)

args = parser.parse_args()

infile = args.input
efile = args.enrich
outfile = args.outfile


# Import files for subsequent analysis
def import_file(filename) -> pd.DataFrame:
    """This function will import and format replicate files."""
    df = pd.read_excel(
        filename,
        # sheet_name="Results",
        # skiprows=46,
        usecols=["Sample Name", "Target Name", "CT"],
        na_values="Undetermined",
    )
    return df


# Import input and enrich files and append any replicate files
files_input = []
for file in infile:
    file_in = import_file(file)
    files_input.append(file_in)
files_enrich = []
for file in efile:
    file_in = import_file(file)
    files_enrich.append(file_in)

# Concat input and enrich into dataframes, then subtract enrich from input
concat_input = pd.concat(files_input, ignore_index=True)
concat_enrich = pd.concat(files_enrich, ignore_index=True)

# Remove _input and _enrich from original data
concat_input["Sample Name"] = concat_input["Sample Name"].str.split(" in").str.get(0)
concat_enrich["Sample Name"] = concat_enrich["Sample Name"].str.split(" en").str.get(0)

ct_recover = pd.DataFrame(
    {
        "Sample Name": concat_input["Sample Name"],
        "Target Name": concat_input["Target Name"],
        "CT_input": concat_input["CT"],
        "CT_enrich": concat_enrich["CT"],
        "CT_diff": concat_input["CT"].sub(concat_enrich["CT"]),
    }
)
ct_recover["Name"] = ct_recover["Sample Name"].str.cat(
    ct_recover["Target Name"], sep=" "
)
ct_recover["percent_recovery"] = 100 * 2 ** ct_recover["CT_diff"]
# Save final table to excel
ct_recover.to_excel(f"{outfile}_output.xlsx")

# Create a barplot and format with all data
# Subsetting data into individual plots can be done with a new script
# using the output file created above
plot = sns.barplot(
    data=ct_recover, x="Name", y="percent_recovery", errorbar="sd", palette="Blues"
)
plt.xticks(rotation=45, ha="right")
plt.subplots_adjust(bottom=0.2)
plt.savefig(f"{outfile}_output.png", dpi=300)
