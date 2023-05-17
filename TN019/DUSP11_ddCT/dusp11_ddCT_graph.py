import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


infile = pd.read_excel('ddCT_DUSP11_table.xlsx')

sns.set(style="ticks")
plot = sns.barplot(data=infile, x='Sample Name', y='value', hue='foldchange', errorbar="sd", palette="Paired")
plt.xlabel(' ')
plt.ylabel('foldchange')
plt.legend(title = '')
plt.savefig('DUSP11_ddCT.png', dpi=300, format="png")