import sys
import pandas as pd


# import 3 replicate files to use later

def replace_char(s, target_char, replacement_char):
    """Replace a character with another character."""
    return s.replace(target_char, replacement_char)


def import_file(filename) -> pd.DataFrame:
    """This function will import, take the mean, and format replicate files."""
    df = pd.read_excel(filename,
                       sheet_name='Results',
                       skiprows=46,
                       usecols=['Sample Name', 'Target Name', 'CT'],
                       na_values = 'Undetermined'
                       )\
        .groupby(by=['Sample Name', 'Target Name'])\
        .mean()\
        .reset_index()
    df['Sample Name'] = df['Sample Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    df['Target Name'] = df['Target Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    return (df)


rep1 = import_file(sys.argv[1])
rep2 = import_file(sys.argv[2])
rep3 = import_file(sys.argv[3])

# Now drop DUSP11KO samples so we only have WT
def drop_row(df: pd.DataFrame, x) -> pd.DataFrame:
    """Drop a row if it contains a string"""
    str_drop = str(x)
    mask = df['Sample Name'].str.contains(str_drop)
    df = df.drop(df[mask].index)
    return (df)

rep1 = drop_row(rep1, 'DUSP11KO')
rep2 = drop_row(rep2, 'DUSP11KO')
rep3 = drop_row(rep3, 'DUSP11KO')

frames = [rep1, rep2, rep3]
combo = pd.concat(frames)

# Get T16 and T24 separated into 18s and DUSP11

t16_18s_subset = combo.loc[combo['Sample Name'].str.contains('WT T16') &
                    combo['Target Name'].str.contains('18S')
]
t16_DUSP11_subset = combo.loc[combo['Sample Name'].str.contains('WT T16') &
                    combo['Target Name'].str.contains('DUSP11')
]
t16_18s_subset.index = range(len(t16_18s_subset))
t16_DUSP11_subset.index = range(len(t16_DUSP11_subset))
# t16_18s_subset= t16_18s_subset[~t16_18s_subset.index.duplicated(keep='first')]
t16 = pd.concat([t16_18s_subset, t16_DUSP11_subset], ignore_index=True, axis=1)
t16.groupby(['0'].str.get(0))
print(t16)