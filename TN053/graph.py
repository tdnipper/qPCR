import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def filter_data(data, filter_value):
    data_new = data[~data['Target Name'].str.contains(filter_value)]
    return data_new

data = filter_data(pd.read_excel('TN053.1_ddCT_data.xlsx'), 'RNA18S1')

def plot_data(data, outfile):
    # Define a custom order for 'Sample Name'
    custom_order = ['0uM', '5uM', '10uM', '25uM']
    data['Sample Name'] = pd.Categorical(data['Sample Name'], categories=custom_order, ordered=True)

    # Sort the data by the custom order
    data = data.sort_values(by='Sample Name')

    # Define unique categories for x-axis and hue
    sample_names = data['Sample Name'].unique()
    target_names = data['Target Name'].unique()

    # Define bar width and positions
    bar_width = 0.35
    x = np.arange(len(sample_names))  # x-coordinates for the groups

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot bars for each target name (hue)
    for i, target in enumerate(target_names):
        subset = data[data['Target Name'] == target]
        ax.bar(
            x + i * bar_width,  # Offset each group by bar_width
            subset['foldchange_mean'],
            bar_width,
            label=target,
            yerr=subset['foldchange_se'],  # Add error bars
            capsize=5
        )

    # Customize the plot
    ax.set_xticks(x + bar_width / 2)  # Center the ticks
    ax.set_xticklabels(sample_names, rotation=45)
    ax.set_ylabel('Fold Change')
    ax.set_title('TN053.1 ddCT Data')
    ax.legend(title='Target Name')

    # Save the plot
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_data(data, 'TN053.1_ddCT_plot_fixed.png')