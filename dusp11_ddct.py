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
                       usecols=['Sample Name', 'Target Name', 'CT']
                       )\
        .groupby(by=['Sample Name', 'Target Name'])\
        .mean()\
        .reset_index()
    df['Sample Name'] = df['Sample Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    df['Target Name'] = df['Target Name'].apply(lambda x: replace_char(x, target_char='_', replacement_char=' '))
    df['Sample Name'] = df['Sample Name'].apply(lambda x: )
    return (df)


rep1 = import_file(sys.argv[1])
rep2 = import_file(sys.argv[2])
rep3 = import_file(sys.argv[3])

def drop_row(df: pd.DataFrame) -> pd.DataFrame:
    """Drop a row if it contains DUSP11KO"""
    str_drop = 'DUSP11KO'
    mask =df['Sample Name'].str.contains(str_drop)
    df = df.drop(df[mask].index)
    return(df)

print(rep2)
# frames = [rep1, rep2, rep3]
# combo = pd.concat(frames)
# print(combo)
