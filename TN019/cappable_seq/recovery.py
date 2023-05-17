import pandas as pd
import os

cwd = os.getcwd()
os.chdir(cwd)


def import_file(filename) -> pd.DataFrame:
    """This function will import and format replicate files."""
    df = pd.read_excel(
        filename,
        sheet_name="Results",
        skiprows=46,
        usecols=["Sample Name", "Target Name", "CT"],
        na_values="Undetermined",
    )
    return df


# Import input and enrich files and append any replicate files
files_input = [
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_11_qPCR1_input.xls",
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_17_t16_input.xls",
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_20_t16_input.xls",
]
input_dfs = []
for file in files_input:
    file_in = import_file(file)
    input_dfs.append(file_in)
files_enrich = [
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_11_enriched.xlsx",
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_17_t16_enriched.xls",
    "/Users/tnipper/Documents/python/qPCR/TN019/cappable_seq/TN019_20_t16_enriched.xls",
]
enrich_dfs = []
for file in files_enrich:
    file_in = import_file(file)
    enrich_dfs.append(file_in)

# Concat input and enrich into dataframes, then subtract enrich from input
concat_input = pd.concat(input_dfs, ignore_index=True)
concat_enrich = pd.concat(enrich_dfs, ignore_index=True)
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
# ct_recover = ct_recover.groupby(by=["Sample Name", "Target Name"])[
#     "percent_recovery"
# ].mean()
# Save final table to excel
ct_recover.to_excel(f"cappable_seq_test_output.xlsx")
