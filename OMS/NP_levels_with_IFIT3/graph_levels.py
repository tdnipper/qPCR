import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

data = pd.read_excel('NP_IFIT3_levels_ddCT.xlsx')
data_filtered = data[data['Target Name'] == 'WSN_NP']

def plot_data(data, outname):
    # Assuming 'foldchange_se' contains the standard error values
    lower_error = data['foldchange_se']
    upper_error = data['foldchange_se']
    
    ax = sns.barplot(
        x='Sample Name', 
        y='foldchange_mean', 
        data=data
    )
    ax.errorbar(
        x=data['Sample Name'], 
        y=data['foldchange_mean'], 
        yerr=[lower_error, upper_error], 
        fmt='none', 
        capsize=5, 
        color='black'
    )
    ax.set_title('NP levels +- IFIT3')
    ax.set_ylabel('Fold Change')
    ax.set_xlabel('')
    custom_labels = ['NP only 1', 'NP only 2', 'NP only 3', 'NP + IFIT3 1', 'NP + IFIT3 2', 'NP + IFIT3 3']
    ax.set_xticklabels(custom_labels, rotation=45)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_data(data_filtered, 'NP_levels_with_IFIT3_fixed.png')