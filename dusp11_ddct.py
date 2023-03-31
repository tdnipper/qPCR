import argparse
import sys
import pandas as pd
import numpy as np
import seaborn as sns

# create args to define variables
parser = argparse.ArgumentParser(description='Analyze ddCT of qPCR files. Can take replicates if listed under inputs flag')

parser.add_argument('-i', '--input', nargs='+', type=str, help='List input qPCR files')
parser.add_argument('-c', '--control', type=str, help='Control or housekeeping gene')
parser.add_argument('-t', '--target', type=str, help='Target gene')
parser.add_argument('-o', '--outfile', type=str, help='File to save ddCT table as')

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
    df = pd.read_excel(filename,
                       sheet_name='Results',
                       skiprows=46,
                       usecols=['Sample Name', 'Target Name', 'CT'],
                       na_values='Undetermined',
                       )\
        .groupby(by=['Sample Name', 'Target Name'])\
        .mean()\
        .reset_index()
    df['Sample Name'] = df['Sample Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    df['Target Name'] = df['Target Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    return (df)

# Use a loop to iterate import_file over all filenames in
files=[]
for file in infile:
    file_in = import_file(file)
    files.append(file_in)
# rep1 = import_file(sys.argv[1])
# rep2 = import_file(sys.argv[2])
# rep3 = import_file(sys.argv[3])
concat = pd.concat(files, ignore_index=True)

# Now drop DUSP11KO samples so we only have WT
def drop_row(df: pd.DataFrame, x) -> pd.DataFrame:
    """Drop a row if it contains a string"""
    str_drop = str(x)
    mask = df['Sample Name'].str.contains(str_drop)
    df = df.drop(df[mask].index)
    return (df)

combo = drop_row(concat, 'DUSP11KO')
print(combo)
# rep1 = drop_row(rep1, 'DUSP11KO')
# rep2 = drop_row(rep2, 'DUSP11KO')
# rep3 = drop_row(rep3, 'DUSP11KO')

# frames = [rep1, rep2, rep3]
# combo = pd.concat(frames)

# Get T16 and T24 separated into 18s and DUSP11

t16_control_subset = combo.loc[combo['Sample Name'].str.contains('WT T16') &
                               combo['Target Name'].str.contains(control)
                               ]
t16_target_subset = combo.loc[combo['Sample Name'].str.contains('WT T16') &
                    combo['Target Name'].str.contains(target)
]

# Combine 18s and DUSP11 data, index according to mock or ix, and take delta CT and mean of each
t16_control_subset.index = range(len(t16_control_subset))
t16_control_subset = t16_control_subset.add_suffix(f'.{control}')
t16_target_subset.index = range(len(t16_target_subset))
t16_target_subset = t16_target_subset.add_suffix(f'.{target}')

t16 = pd.concat([t16_control_subset, t16_target_subset],
                axis=1,
                )
# Use combined data to get deltaCT and clean up dataframe
t16['delta'] = t16[f'CT.{target}']-t16[f'CT.{control}']
t16 = t16.drop([f'Target Name.{control}', f'CT.{control}', f'Sample Name.{target}', f'Target Name.{target}', f'CT.{target}'], axis=1)

# Use Sample Name column as index to group by mock and ix to get mean and std of deltaCT
t16 = t16.rename(columns={f'Sample Name.{control}': 'Sample Name'})
t16.set_index('Sample Name', inplace=True)
t16['mean'] = t16.groupby(by=['Sample Name'])['delta'].mean()
t16['sd'] = t16.groupby(by=['Sample Name'])['delta'].std()
t16.reset_index(inplace=True)

# Separate out mock and ix, then subtract ix-mock to get deltadeltaCT
t16_mock = t16.loc[t16['Sample Name'] == 'WT T16 mock']\
    .drop(['delta'], axis=1)\
    .drop_duplicates()\
    .add_suffix('.mock')\
    .reset_index()
t16_ix = t16.loc[t16['Sample Name'] == 'WT T16 ix']\
    .drop(['delta'], axis=1)\
    .drop_duplicates()\
    .add_suffix('.ix')\
    .reset_index()
t16_final = pd.concat([t16_mock, t16_ix], axis=1)
t16_final['ddCT'] = t16_final['mean.mock']-t16_final['mean.ix']
t16_final['foldchange'] = np.power(2, t16_final['ddCT'])
print(t16_final)