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
        # usecols=["Sample Name", "Target Name", "CT"],
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

# Group based on MOI and timepoint
ct_recover["group"] = (
    ct_recover["Sample Name"]
    .astype(str)
    .apply(lambda x: "MOI5 T10" if x.startswith("MOI5 T10") else "MOI0.02 T24")
)

# Assign infection status based on mock or ix in sample name
ct_recover["infect"] = (
    ct_recover["Sample Name"]
    .astype(str)
    .apply(lambda x: "mock" if x.endswith("mock") else ("ix"))
)

# Save final table to excel
ct_recover.to_excel(f"{outfile}_output.xlsx")

colors = ["#00a9e0", "#00549f"]
# Create a barplot and format with all data
# Subsetting data into individual plots can be done with a new script
# using the output file created above
g = sns.catplot(
    data=ct_recover,
    kind="bar",
    x="Target Name",
    y="percent_recovery",
    hue="infect",
    col="group",
    palette=colors,
    col_wrap=4,
    height=4,
    aspect=1,
    legend=False,
)

# Access the Axes objects for each subplot
axes = g.axes.flatten()

# Rotate the x-axis labels by 45 degrees on both subplots
for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
# plt.xticks(rotation=45, ha="right")
plt.subplots_adjust(left=0.2, bottom=0.4)
plt.savefig(f"{outfile}_output.png", dpi=300)
plt.savefig(f"{outfile}_output.svg", dpi=300)

g.set(ylim=(0, 15))
plt.savefig(f"{outfile}_output_zoom.svg", dpi=300)
plt.savefig(f"{outfile}_output_zoom.png", dpi=300)
