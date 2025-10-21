import pandas as pd

data = pd.read_excel('TN061_3.xlsx', sheet_name='Results', skiprows=46, usecols=['Sample Name', 'Target Name', 'CT'], na_values=['Undetermined', 'NTC']).reset_index(drop=True).dropna()
data.to_excel('filtered_results.xlsx', index=False)