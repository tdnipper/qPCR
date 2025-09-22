import pandas as pd

data = pd.read_excel('TN052_13_redo.xlsx', sheet_name='Results', skiprows=46, usecols=['Sample Name', 'Target Name', 'Task', 'CT'], na_values=['Undetermined', 'NTC']).reset_index(drop=True).dropna()
data['Task'] = data['Task'].astype(str)
# data['Sample Name'] = data['Sample Name'] + '_' + data['Task']
data.to_excel('filtered_results.xlsx', index=False)