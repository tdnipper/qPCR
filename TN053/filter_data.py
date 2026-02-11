import pandas as pd

def filter_data(input_file, output_file):
    # Read the CSV file
    df = pd.read_excel(input_file)

    # Filter the DataFrame to include only rows where 'column_name' is greater than 0
    filtered_df = df[(~df['Sample Name'].str.contains('noRT')) & (df['Sample Name'] != 'NTC')]

    # Save the filtered DataFrame to a new CSV file
    filtered_df.to_excel(output_file, index=False)

if __name__ == "__main__":
    input_file = 'TN053.1_filtered.xlsx'
    output_file = 'TN053.1_fixed.xlsx'
    filter_data(input_file, output_file)